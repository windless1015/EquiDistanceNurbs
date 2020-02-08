#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_EquiDistanceNurbs.h"
#include <qvector.h>

//最大值, 最小值
const double DMAX = 1.e12;
const double DMIN = 1.e-12;

typedef struct WeightPoint
{
	QPointF point;
	double weight;
} structPoint;


class EquiDistanceNurbs : public QMainWindow
{
	Q_OBJECT

public:
	EquiDistanceNurbs(QWidget *parent = Q_NULLPTR);



public slots:
	void OnBtnEquidistance();

protected:
	void mousePressEvent(QMouseEvent* event);
	void paintEvent(QPaintEvent* event);

private:
	Ui::EquiDistanceNurbsClass ui;
	//三阶 样条
	const int DEGREE = 3;
	//点之间的步长距离
	const double step = 0.01;

	//原始nurbs曲线的控制点
	QVector<structPoint>			nurbsCtrlPoints;
	//原始nurbs曲线上的点
	QVector<QPointF>			    originNurbsBodyPts;

private:
	bool isShowCtrlPoints; //是否显示控制点
	bool isShowCtrlPtsConnectedLine; //是否显示控制点之间的直连线

private:
	void GenerateNURBSCurve();
	double N3Spline(int i, double u);

};
