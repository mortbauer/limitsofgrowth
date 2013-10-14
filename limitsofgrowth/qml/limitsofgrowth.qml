import QtQuick 2.0

Item {
    id:container
    width: 320
    height: 480

    Column {
        spacing:5
        anchors.fill: parent
        anchors.topMargin: 12

        Text {
            font.pointSize: 12
            font.bold: true
            text: "Limits of Growth, simple World Model"
            anchors.horizontalCenter: parent.horizontalCenter
            color: "#777"
        }

        Canvas {
            id:canvas
            width:320
            height:280
            property color strokeStyle:  Qt.darker(fillStyle, 1.4)
            property color fillStyle: "#b40000" // red
            property int lineWidth: lineWidthCtrl.value
            property bool fill: true
            property bool stroke: true
            property real alpha: 1.0
            property real scale : scaleCtrl.value
            property real rotate : rotateCtrl.value
            antialiasing: true

            onLineWidthChanged:requestPaint();
            onScaleChanged:requestPaint();
            onRotateChanged:requestPaint();

            onPaint: {
                var ctx = canvas.getContext('2d');
                var originX = 85
                var originY = 75
                ctx.save();
                ctx.clearRect(0, 0, canvas.width, canvas.height);
                ctx.translate(originX, originX);
                ctx.globalAlpha = canvas.alpha;
                ctx.strokeStyle = canvas.strokeStyle;
                ctx.fillStyle = canvas.fillStyle;
                ctx.lineWidth = canvas.lineWidth;

                ctx.translate(originX, originY)
                ctx.scale(canvas.scale, canvas.scale);
                ctx.rotate(canvas.rotate);
                ctx.translate(-originX, -originY)

                //! [0]
                ctx.beginPath();
                for (var i = 0; i < 100; i +=1){
                    ctx.moveTo(i,world.population[i]);
                }
                ctx.closePath();
                //! [0]
                if (canvas.fill)
                    ctx.fill();
                if (canvas.stroke)
                    ctx.stroke();
                ctx.restore();
            }
        }
    }
    Column {
        id: controls
        anchors.bottom: parent.bottom
        anchors.bottomMargin: 12
        Slider {id: lineWidthCtrl; min: 1; max: 10; init: 2; name: "Outline"}
        Slider {id: scaleCtrl; min: 0.1; max: 10; init: 1; name: "Scale"}
        Slider {id: rotateCtrl; min: 0; max: Math.PI*2; init: 0; name: "Rotate"}
    }
}
