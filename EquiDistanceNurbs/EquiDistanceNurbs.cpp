#include "EquiDistanceNurbs.h"
#include <qevent.h>
#include <qpainter.h>
#include <qtextstream.h>
#include <qdebug.h>

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
		NurbsCtrlPoint oneCtrlPt;
		oneCtrlPt.point = event->pos();;
		oneCtrlPt.weight = 1.0;
		nurbsCtrlPoints.push_back(oneCtrlPt);
		GenerateNURBSCurve();
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
		ctrlPtsPen.setWidth(5);
		painter.setPen(ctrlPtsPen);
		QBrush brush(Qt::SolidPattern);//画刷
		brush.setColor(Qt::blue);
		painter.setBrush(brush);//设置画刷
		for (int i = 0; i < nurbsCtrlPoints.size(); i++)
		{
			//绘制点使用drawEllipse ,把ctrlPtsPen.setWidth设置为5(或更大) 点的显示更明显
			painter.drawEllipse(static_cast<int>(nurbsCtrlPoints[i].point.x()) ,
				static_cast<int>(nurbsCtrlPoints[i].point.y()) ,2, 2);
			//painter.drawPoint(nurbsCtrlPoints[i].point.x(), nurbsCtrlPoints[i].point.y());
		}
	}

	//显示控制点之间直连线
	if (isShowCtrlPtsConnectedLine)
	{
		QPen ctrlPen(QColor(0, 255, 0));
		ctrlPen.setWidth(2);
		ctrlPen.setStyle(Qt::DotLine); //设置直线样式为点虚线
		painter.setPen(ctrlPen);
		for (int i = 0; i < nurbsCtrlPoints.size() - 1; i++)
		{
			painter.drawLine(nurbsCtrlPoints[i].point, nurbsCtrlPoints[i + 1].point);
		}
	}

	//绘制NURBS曲线
	QPen curvePen(QColor(0, 0, 0));
	curvePen.setWidth(3);
	painter.setPen(curvePen);

	for (int i = 0; i < originNurbsBodyPts.size() - 1; i++) //遍历body pts是从idx为[0, i-2]
	{
		//绘制每一小段的直线,当body pts足够多的时候就是曲线
		painter.drawLine(originNurbsBodyPts[i].point, originNurbsBodyPts[i + 1].point);
	}

}

double EquiDistanceNurbs::N3Spline(int i, double u)
{
	double t = u - i;
	double res;
	if (0 <= t && t < 1)
		res =  t*t*t;
	else if (1 <= t && t < 2) 
		res = (-3 * pow(t - 1, 3) + 3 * pow(t - 1, 2) + 3 * (t - 1) + 1);
	else if (2 <= t && t < 3)
		res = (3 * pow(t - 2, 3) - 6 * pow(t - 2, 2) + 4);
	else if (3 <= t && t < 4)
		res = pow(4 - t, 3);
	else
		res = 0;
	return 1.0 / 6.0 * res;
}


void EquiDistanceNurbs::GenerateNURBSCurve()
{
	originNurbsBodyPts.clear();
	//int count = 0;
	for (double u = 1; u < nurbsCtrlPoints.size() + 2; u += step)
	{
		NurbsBodyPoint bodyPt;
		double factor = 0;
		double N;
		for (int i = 0; i < nurbsCtrlPoints.size(); i++)
		{
			N = N3Spline(i, u);
			bodyPt.point += nurbsCtrlPoints[i].point * N;
			factor += N;
		}
		bodyPt.point /= factor;
		originNurbsBodyPts.push_back(bodyPt);
		//count++;
	}
	//qDebug() << "count : " << count << endl;
}

