import logging
import pandas
import numpy as np
from scipy.integrate import ode, odeint
from scipy import optimize
from numpy import linalg
from collections import OrderedDict
import glob
import json
import traversedata
import requests
from math import log
import re
from pkg_resources import resource_filename

logger = logging.getLogger(__name__)

class WorldSimple(object) :
    params = OrderedDict((
        ('birthrate',{'initial':0.03,'max':10,'min':0.01}),
        ('deathrate',{'initial':0.01,'max':10,'min':0.01}),
        ('regenerationrate',{'initial':0.1,'max':10,'min':0.01}),
        ('burdenrate',{'initial':0.02,'max':10,'min':0.01}),
        ('economyaim',{'initial':10,'max':10,'min':0.1}),
        ('growthrate',{'initial':0.05,'max':10,'min':0.01})
    ))

    cols = ['population','burden','economy']

    reference_data_names = {'population':'SP.POP.TOTL',
                            'economy':'NY.GDP.MKTP.CD',
                            'burden':'EN.ATM.CO2E.KT',
                            }

    reference_parameters = {'birthrate':{'name':'SP.DYN.CBRT.IN','per':1000},
                            'deathrate':{'name':'SP.DYN.CDRT.IN','per':1000},
                            'burdenrate':{'name':'EN.ATM.CO2E.KT','per':'SP.POP.TOTL'},
                            }

    time = np.arange(1900,1900+250+1,1)

    def __init__(self, resolution=250):
        self.sensitivities = []
        self.initial = {param:v['initial'] for param,v in self.params.items()}
        for column in self.cols:
            for param in self.params:
                self.sensitivities.append('{0},{1}'.format(column,param))

    def get_ref_parameters(self,country,year=1960):
        d = {}
        for param,v in self.reference_parameters.items():
            if isinstance(v['per'],str):
                per = self._reference_data(country)[v['per']]
            else:
                per = v['per']
            d[param] = (self._reference_data(country)[v['name']]/per).loc[year]
        return d


    def create_input_vector(self,params):
        return [params[param] for param in self.params]

    @staticmethod
    def birth_mod(p,x):
        fac = 1-(x[2]-1)/p[4]
        if fac >=0:
            return p[0]*fac*x[0]/x[1]*x[2]
        else:
            return p[0]*0.1*x[0]/x[1]*x[2]

    @staticmethod
    def create_dx(p,birth_callback=None):
        """parameter vector p must be in correct order, to ensure this you can
        construct it throught the function create_input_vector which takes a
        dictionary as input
        """
        birthrate,deathrate,regenerationrate,burdenrate,economyaim,growthrate=p
        def func(time,x):
            """
            the change of the system variables population, burden and economy
            x: [population,burden,economy]
            """
            population,burden,economy = x
            quality = burden**(-1)
            if not birth_callback:
                birth = birthrate * population * quality * economy
            else:
                birth = birth_callback(p=p,x=x)
            death = population * deathrate * burden
            ecocide = burdenrate * economy * population
            if quality > 1:
                regeneration = regenerationrate * burden
            else:
                regeneration = regenerationrate
            economicgrowth = growthrate * economy * burden \
                    * (1-(economy*burden)/economyaim)
            f = lambda : [birth-death,ecocide-regeneration,economicgrowth]
            fnachx = lambda : np.array([
                [birthrate/burden*economy-deathrate*burden,
                birthrate*population*(-1)/burden**2*economy-population*deathrate,
                birthrate*population/burden],
                [burdenrate*economy,
                regenerationrate if burden < 1 else 0,
                burdenrate*population],
                [0,
                growthrate*economy-2*growthrate*economy**2*burden/economyaim,
                growthrate*burden-2*growthrate*burden**2*economy/economyaim]
            ])
            fnachp = lambda : np.array([
                [population/burden*economy,
                - population+burden,0,0,0,0],
                [0,0,-burden if burden<1 else -1,
                economy*population,0,0],
                [0,0,0,0,growthrate*economy**2*burden**2/economyaim**2,
                economyaim*burden*(1-(economy*burden/economyaim))]
            ])
            return {'f':f,'fx':fnachx,'fp':fnachp}
        return func

    @staticmethod
    def create_dx_mod(p_hic,p_lic):
        """parameter vector p must be in correct order, to ensure this you can
        construct it throught the function create_input_vector which takes a
        dictionary as input
        """
        birthrate_hic,deathrate_hic,regenerationrate_hic,burdenrate_hic,economyaim_hic,growthrate_hic=p_hic
        birthrate_lic,deathrate_lic,regenerationrate_lic,burdenrate_lic,economyaim_lic,growthrate_lic=p_lic
        def func(time,x):
            """
            the change of the system variables population, burden and economy
            x: [population,burden,economy]
            """
            population_hic,burden_hic,economy_hic,population_lic,burden_lic,economy_lic = x
            quality_hic = burden_hic**(-1)
            quality_lic = burden_lic**(-1)
            birth_hic = birthrate_hic * population_hic * quality_hic * economy_hic
            birth_lic = birthrate_lic * population_lic * quality_lic * economy_lic
            death_hic = population_hic * deathrate_hic * burden_hic
            death_lic = population_lic * deathrate_lic * burden_lic
            ecocide_hic = burdenrate_hic * economy_hic * population_hic
            ecocide_lic = burdenrate_lic * economy_lic * economy_hic * population_lic
            if quality_hic > 1:
                regeneration_hic = regenerationrate_hic * burden_hic
            else:
                regeneration_hic = regenerationrate_hic
            if quality_lic > 1:
                regeneration_lic = regenerationrate_lic * burden_lic
            else:
                regeneration_lic = regenerationrate_lic
            economicgrowth_hic = growthrate_hic * economy_hic * burden_hic \
                    * (1-(economy_hic*burden_hic)/economyaim_hic)
            economicgrowth_lic = growthrate_lic * economy_lic * burden_lic \
                    * (1-(economy_lic*burden_lic)/economyaim_lic)
            return [birth_hic-death_hic,ecocide_hic-regeneration_hic,economicgrowth_hic,
                    birth_lic-death_lic,ecocide_lic-regeneration_lic,economicgrowth_lic]
        return func

    def x_odeint(self,params={},func=None):
        if not func:
            func = self.create_dx(self.create_input_vector(params))
        res,info = odeint(lambda x,t:func(t,x)['f'](),[1.0,1.0,1.0], self.time,
                          full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            dataframe = pandas.DataFrame(res,
                columns=['population','burden','economy'],index=self.time)
            return dataframe

    def x_odeint_mod(self,params=(),func=None,x0={}):
        if not func:
            func = self.create_dx_mod(self.create_input_vector(params[0]),
                                      self.create_input_vector(params[1]))
        if not x0:
            x0 = np.ones(6)
        res,info = odeint(lambda x,t:func(t,x),x0, self.time,
                          full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            dataframe_hic = pandas.DataFrame(res[:,:3],
                columns=['population','burden','economy'],index=self.time)
            dataframe_lic = pandas.DataFrame(res[:,3:],
                columns=['population','burden','economy'],index=self.time)
            return dataframe_hic,dataframe_lic

    def x_cvode(self,params):
        from assimulo.problem import Explicit_Problem
        from assimulo.solvers import CVode
        func = self.create_dx(self.create_input_vector(params))
        problem = Explicit_Problem(lambda t,x:func(t,x)['f'](), [1.0,1.0,1.0],0)
        sim = CVode(problem)
        t,x = sim.simulate(250,len(self.time)-1)
        dataframe = pandas.DataFrame(x,
                columns=['population','burden','economy'],index=self.time)
        return dataframe

    def create_ds(self,p):
        f = self.create_dx(p)
        dx = np.empty((21,))
        def func(t,x):
            _f = f(t,x[:3])
            dx[:3] = _f['f']()
            _s = x[3:].reshape((3,6))
            dx[3:] = (_f['fx']().dot(_s)+_f['fp']()).reshape((18,))
            return dx
        return func

    def s_odeint(self,params):
        func = self.create_ds(self.create_input_vector(params))
        s0 = np.ones(21)
        s0[3:]=0
        res,info = odeint(lambda s,t:func(t,s),s0, self.time,
                          full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            dataframe = pandas.DataFrame(
                res,columns=self.cols+self.sensitivities,index=self.time)
            return dataframe

    def s_cvode(self,params):
        from assimulo.problem import Explicit_Problem
        from assimulo.solvers import CVode
        func = self.create_ds(self.create_input_vector(params))
        s0 = np.ones(21)
        s0[3:]=0
        problem = Explicit_Problem(lambda t,s:func(t,s),s0,1900)
        sim = CVode(problem)
        t,s = sim.simulate(1900+250,self.time.shape[0]-1)
        dataframe = pandas.DataFrame(
            s,columns=self.cols+self.sensitivities,index=self.time)
        return dataframe

    def s_cvode_natural(self,params):
        from assimulo.problem import Explicit_Problem
        from assimulo.solvers import CVode
        problem = Explicit_Problem(lambda t,x,p:self.create_dx(p)(t,x)['f'](),
                                   [1.0,1.0,1.0],0,
                                   [params[p] for p in self.params])
        sim = CVode(problem)
        sim.report_continuously = True
        t,x = sim.simulate(250,self.time.shape[0]-1)
        dataframe = pandas.DataFrame(x,
                columns=['population','burden','economy'])
        d = {}
        sens = np.array(sim.p_sol)
        for i,col in enumerate(self.cols):
            for j,param in enumerate(
                ('birthrate','deathrate','regenerationrate',
                 'burdenrate','economyaim','growthrate')):
                d['{0},{1}'.format(col,param)] = sens[j,:,i]
        dataframe_sens = pandas.DataFrame(d,index=self.time)
        return dataframe_sens

    def s_ode(self,params):
        x0 = np.ones(3)
        s0 = np.zeros(18)
        solver_x = ode(self.dx).set_integrator('dopri5')
        solver_x.set_initial_value(x0,0).set_f_params(params,)

        solver_s = ode(self.ds).set_integrator('dopri5')
        solver_s.set_initial_value(s0,0).set_f_params(params,x0)

        dt = 1
        t1 = 250
        sensitivities = []
        for column in self.cols:
            for param in self.params:
                sensitivities.append('{0},{1}'.format(column,param))
        sol = pandas.DataFrame(np.empty((t1/dt,22)),columns=self.cols+sensitivities)
        i = 0
        #return solver_x,solver_s
        while solver_x.successful() and solver_s.successful() and solver_x.t < t1:
            solver_x.integrate(solver_x.t+dt)
            sol.iloc[i][self.cols] = solver_x.y
            solver_s.set_f_params(params,solver_x.y)
            solver_s.integrate(solver_s.t+dt)
            sol.iloc[i][sensitivities] = solver_s.y
            i += 1
        return sol

    def _reference_data(self,country):
        if not hasattr(self,'_wc'):
            wc = WorldBankClient()
            wc.load_from_hdf(resource_filename(__name__,'data/world_indicators.hdf'))
            self._wc = wc
        d = self._wc.indicators_by_countries(country)
        d.index = d.index.astype(np.int)
        return d

    def reference_data(self,country):
        return self._reference_data(country).rename_axis(
            {v:key for key,v in self.reference_data_names.items()})


    def create_bnds(self):
        return [(self.params[p]['min'],self.params[p]['max']) for p in self.params]

    def create_residum_func(self,ref,weights={
        'economy':10,'population':1,'burden':1},x_callback=None):
        n = self.reference_data('WLD').shape[0]
        d = np.empty(n*3)
        if not x_callback:
            res = lambda p:self.x_odeint(func=self.create_dx(p))
        else:
            res = lambda p:x_callback(p)
        norm_ref = ref/ref.loc[1960]
        def residum(p):
            try:
                r = res(p)
                diff = r.loc[norm_ref.index]/r.loc[1960]-norm_ref
                for i,c in enumerate(diff):
                    d[i*n:i*n+n] = (diff[c]*weights[c])**2
            except Exception as e:
                d[:] += 10e6
                print('######',e)
                #logger.warn('failed to calc diff with p = {0}'.format(p))
            return d
        return residum

    def fit_with_data_basinhopping(self,country='WLD'):
        func = self.create_residum_func(self.reference_data(country))
        _min = np.array([bnd[0] for bnd in bnds])
        _max = np.array([bnd[1] for bnd in bnds])

        def accept_test(**kwargs):
            x = kwargs['x_new']
            res =  bool(np.all(x<=_max)) and bool(np.all(x>=_min))
            return res

        minimizer_kwargs = {"method":"L-BFGS-B", "jac":False,'bounds':self.create_bnds()}
        r = optimize.basinhopping(
            lambda *args:linalg.norm(func(*args)),
            self.create_input_vector(self.initial),accept_test=accept_test,
            niter=100,minimizer_kwargs=minimizer_kwargs
        )

        return r

    def fit_with_data_bfgs(self,country='WLD',**kwargs):
        func = self.create_residum_func(self.reference_data(country),**kwargs)
        r = optimize.fmin_l_bfgs_b(
            lambda *args:linalg.norm(func(*args)),self.create_input_vector(self.initial),
            bounds=self.create_bnds(),approx_grad=True,pgtol=1e-6
        )
        return r

    def fit_with_data_bfgs_mod(self,country='WLD',**kwargs):
        w = self
        rel = w.reference_data('HIC').loc[1960]/w.reference_data('LMY').loc[1960]
        x0_lic = (rel+1)**(-1)
        x0_hic = x0_lic*rel
        p_hic = w.initial.copy()
        p_lic = w.initial.copy()
        x0 = [x0_hic[c] for c in w.cols]+[x0_lic[c] for c in w.cols]

        def callback(p):
            x0[2]=p[0]
            x0[4]=p[1]
            x0[5]=p[3]
            x0[8]=p[4]
            x0[10]=p[5]
            x0[11]=p[6]
            f = self.create_dx_mod(x0[:6],x0[6:])
            x_hic,x_lic = self.x_odeint_mod(func=f)
            return x_hic+x_lic
        func = self.create_residum_func(self.reference_data(country),x_callback=callback,**kwargs)
        r = optimize.fmin_l_bfgs_b(
            lambda p:linalg.norm(func(p)),np.ones(6),
            bounds=self.create_bnds()*2,approx_grad=True,pgtol=1e-6
        )
        return r

    def fit_with_data_leastsq(self,country='WLD'):
        func = self.create_residum_func(self.reference_data(country))
        r = optimize.leastsq(func,self.create_input_vector(self.initial),full_output=True)
        return r

    def co2_contentration(self):
        d = pandas.read_csv(
            resource_filename(__name__,'data/monthly_mlo.csv'),
            skiprows=54,header=[2,0,1],index_col=0,
        )
        return d.loc[self.reference_data('WLD').index].icol(-1).groupby(level=0).mean()
        #return d[d[(' seasonally', 'adjusted filled', '    [ppm]')]>0][(' seasonally', 'adjusted filled', '    [ppm]')]

class WorldBankClient(object):
    _REGION_CODES = 'http://worldbank.270a.info/classification/region.html'
    BASE_URL = 'http://api.worldbank.org/'
    _ID = re.compile(r'world_bank_data_(?P<country>[A-Z]*)_(?P<indicator>[A-Z.0-9]*).json')

    def __init__(self):
        self.indicators = {}
        self._plots = {}

    def indicators_by_countries(self,country):
        indicators = self.indicators.keys()
        d = {ind:v[country] for ind,v in self.indicators.items() if country in v}
        return pandas.DataFrame(d)

    @property
    def all_countries(self):
        if not hasattr(self,'_all_countries'):
            self._all_countries = {x['name']:x['id'] for x in self.countries_basicinfo}
        return self._all_countries

    @property
    def countries_basicinfo(self):
        if not hasattr(self,'_countries_basicinfo'):
            self._countries_basicinfo = self._get(
                '{0}countries/all/'.format(self.BASE_URL))
        return self._countries_basicinfo

    @property
    def indicators_basicinfo(self):
        if not hasattr(self,'_indicators_basicinfo'):
            self._indicators_basicinfo = self._get(
                '{0}indicators/all/'.format(self.BASE_URL))
        return self._indicators_basicinfo

    @property
    def all_indicators(self):
        if not hasattr(self,'_all_indicators'):
            self._all_indicators = {x['name']:x['id'] for x in self.indicators_basicinfo}
        return self._all_indicators

    def get_indicator_by_country(self,indicator,country):
        if not indicator in self.indicators or not country in self.indicators[indicator]:
            d = traversedata.WorldBankData(
                self._get_indicator_by_country(indicator,country),
                country=country,indicator=indicator)
            if indicator in self.indicators:
                self.indicators[indicator][country] = d.data
            else:
                self.indicators[indicator] = pandas.DataFrame(d.data,columns=[country])
        return self.indicators[indicator][country]

    def load_from_json_and_dir(self,path):
        for filepath in glob.glob('{0}world_bank_data_*.json'.format(path)):
            with open(filepath,'r') as f:
                m = self._ID.match(f.name.split('/')[-1])
                indicator = m.group('indicator')
                country = m.group('country')
                if not indicator in self.indicators or not country in self.indicators[indicator]:
                    d = traversedata.WorldBankData(
                        json.loads(f.read()))
                    if indicator in self.indicators:
                        self.indicators[indicator][country] = d
                    else:
                        self.indicators[indicator] = pandas.DataFrame(d.data,columns=[country])

    def save_to_json(self,path):
        with open(path,'w') as f:
            f.write(json.dumps({ind:v.raw for ind,v in self.indicators.items()}))

    def save_to_hdf(self,path,mode='a'):
        store = pandas.HDFStore(path,mode=mode)
        try:
            for indicator,countries in self.indicators.items():
                countries.to_hdf(store,'{0}'.format(
                    indicator.replace('.','_')))
        finally:
            store.close()

    def load_from_hdf(self,path):
        store = pandas.HDFStore(path,mode='r')
        try:
            for key in store.keys():
                indicator = key.lstrip('/')
                indicator = indicator.replace('_','.')
                self.indicators[indicator] = store[key]
        finally:
            store.close()

    def _get_indicator_by_country(self,indicator,country):
        return self._get('{0}countries/{1}/indicators/{2}'.format(
            self.BASE_URL,country,indicator))

    def _get(self,url,page=1,per_page=50,until=None):
        result = []
        def get(page):
            r = requests.get(
                '{0}?format=json&per_page={1}&page={2}'.format(url,per_page,page))
            if r.ok:
                status,res = r.json()
                result.extend(res)
                if page < status['pages'] and until and page < until:
                    get(page+1)
        get(1)
        return result

    def plot(self):
        self._plots['main'] = DataFramePlot()
        self._plots['main']

