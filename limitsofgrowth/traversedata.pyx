import os
import json
import numpy
import pandas

class WorldBankData(object):
    def __init__(self,list data,indicator='',country=''):
        self.raw = data
        self.indicator = indicator
        self.country = country

    @property
    def data(self):
        cdef int i
        if not hasattr(self,'_data'):
            data = numpy.empty(len(self.raw))
            index = numpy.empty(len(self.raw))
            for i in range(len(self.raw)):
                data[i] = self.raw[i]['value']
                index[i] = self.raw[i]['date']
            self._data = pandas.TimeSeries(data,index=index)
        return self._data

    def save_json(self,filename='',path=''):
        if not filename:
            filename = 'world_bank_data_{0}_{1}.json'.format(
                self.country,self.indicator)
        if not path:
            fullpath = os.path.abspath(filename)
        else:
            fullpath = os.path.abspath(os.path.join(path,filename))

        with open(fullpath,'w') as f:
            f.write(json.dumps(self.raw))

    def save_hdf(self,filename='',path='',append=True):
        if not filename:
            filename = 'world_bank_data_{0}_{1}.hdf'.format(
                self.country,self.indicator)
        if not path:
            fullpath = os.path.abspath(filename)
        else:
            fullpath = os.path.abspath(os.path.join(path,filename))

        with pandas.HDFStore(fullpath,mode='a' if append else 'w') as f:
            pass
        # don't know much about hdf yet needs more investigation

