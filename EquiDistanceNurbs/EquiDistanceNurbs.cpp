#include "EquiDistanceNurbs.h"
#include "qevent.h"
#include "qpainter.h"
#include "qtextstream.h"

EquiDistanceNurbs::EquiDistanceNurbs(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	isShowCtrlPoints = true;
	isShowCtrlPtsConnectedLine = true;
	connect(ui.btnOffsetCurve, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnEquidistance);
}


void EquiDistanceNurbs::OnBtnEquidistance()
{

}

void EquiDistanceNurbs::mousePressEvent(QMouseEvent* event)
{
	if (event->button() == Qt::LeftButton)
	{
		structPoint oneCtrlPt;
		oneCtrlPt.point = event->pos();;
		oneCtrlPt.weight = 1.0;
		nurbsCtrlPoints.push_back(oneCtrlPt);
		//generateCurve(Isseal = false);
		update();
	}
}

void EquiDistanceNurbs::paintEvent(QPaintEvent *event)
{
	// draw Contorl Points
	QPainter painter(this);
	if (isShowCtrlPoints)
	{
		if (nurbsCtrlPoints.size() <= 0)
			return;
		QPen ctrlPtsPen(QColor(0, 0, 255));
		ctrlPtsPen.setWidth(7);
		painter.setPen(ctrlPtsPen);
		QBrush brush(Qt::SolidPattern);//画刷
		brush.setColor(Qt::blue);
		painter.setBrush(brush);//设置画刷
		for (int i = 0; i < nurbsCtrlPoints.size(); i++)
		{
			painter.drawEllipse(static_cast<int>(nurbsCtrlPoints[i].point.x()) ,
				static_cast<int>(nurbsCtrlPoints[i].point.y()) ,2, 2);
			//painter.drawPoint(ctrlPoints[i].point.x(), ctrlPoints[i].point.y()); //qt自身接口绘点不好看
		}
	}

	//显示控制点之间直连线
	if (isShowCtrlPtsConnectedLine)
	{
		QPen ctrlPen(QColor(0, 255, 0));
		ctrlPen.setWidth(2);
		ctrlPen.setStyle(Qt::DashDotDotLine);
		painter.setPen(ctrlPen);
		for (int i = 0; i < nurbsCtrlPoints.size() - 1; i++)
		{
			painter.drawLine(nurbsCtrlPoints[i].point, nurbsCtrlPoints[i + 1].point);
		}
	}

	// draw NURBS Curve
	/*QPen curvePen(QColor(0, 0, 0));
	curvePen.setWidth(2);
	painter.setPen(curvePen);

	for (int i = 0; i < curvePoints.size() - 1; i++)
	{
		painter.drawLine(curvePoints[i], curvePoints[i + 1]);
	}*/

	// draw OFFSET Curve
	/*QPen offsetPen(QColor(255, 0, 0));
	offsetPen.setWidth(2);
	painter.setPen(offsetPen);

	for (int i = 0; i < db_offsetPoints.size(); i++)
	{
		for (int j = 0; j < db_offsetPoints.at(i).size() - 1; j++)
		{
			painter.drawLine(db_offsetPoints.at(i)[j], db_offsetPoints.at(i)[j + 1]);
		}
	}

	for (int i = 0; i < offsetPoints.size() - 1; i++)
	{
		painter.drawLine(offsetPoints.at(i), offsetPoints.at(i + 1));
	}*/

}





double EquiDistanceNurbs::N3Spline(int i, double u)
{
	double t = u - i;
	double a = 1.0 / 6.0;

	if (0 <= t && t < 1)
		return a * t*t*t;
	else if (1 <= t && t < 2) 
		return a * (-3 * pow(t - 1, 3) + 3 * pow(t - 1, 2) + 3 * (t - 1) + 1);
	else if (2 <= t && t < 3)
		return a * (3 * pow(t - 2, 3) - 6 * pow(t - 2, 2) + 4);
	else if (3 <= t && t < 4)
		return a * pow(4 - t, 3);
	else return 0;
}


