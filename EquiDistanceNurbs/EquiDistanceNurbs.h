#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_EquiDistanceNurbs.h"
#include <qvector.h>

//定义最小值,最大值
const double DMAX = 1.e12;
const double DMIN = 1.e-12;

typedef struct WeightPoint
{
	QPointF point;
	double weight;
} NurbsCtrlPoint;

typedef struct BodyPoint
{
	QPointF point; //点的坐标
	QPointF firstDerivative; //一阶导数, x表示x方向的一阶导数, y表示y方向的一阶导数
	QPointF secondDerivative; //二阶导数,x表示x方向的二阶导数, y表示y方向的二阶导数
	QPointF unitNormVector; //单位法向量
	double curvatureRadius; //曲率半径
} NurbsBodyPoint;


class EquiDistanceNurbs : public QMainWindow
{
	Q_OBJECT

public:
	EquiDistanceNurbs(QWidget *parent = Q_NULLPTR);



public slots:
	void OnBtnRawEquidistance();
	void OnBtnDisplayTangent();
	void OnBtnDisplayNorm();
	void OnBtnDisplayCurvature();
	void OnDisplayCurvature();


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
	QVector<NurbsCtrlPoint>			nurbsCtrlPoints;
	//原始nurbs曲线上的点
	QVector<NurbsBodyPoint>		originNurbsBodyPts;
	QVector<NurbsBodyPoint> locationVec;

private:
	bool isShowCtrlPoints; //是否显示控制点
	bool isShowCtrlPtsConnectedLine; //是否显示控制点之间的直连线
	bool isShowBodyPtsTangent; //是否显示曲线点的切向量
	bool isShowBodyPtsNorm; //是否显示曲线点的法向量
	bool isShowCurvatureRadius; //是否显示某一个点的曲率半径
	bool isShowRawEquiDistance; //是否显示原始等距线

private:
	void GenerateNURBSCurve();
	double N3Spline(int i, double u);
	//计算一阶导数和法向量
	void CalFirstDerivativeAndNorm(QVector<NurbsBodyPoint>&);
	//计算二阶导数和曲率半径
	void CalSecondDerivativeAndCurvRad(QVector<NurbsBodyPoint>&);
	//计算偏移曲线
	void CalOffsetCurve(const QVector<NurbsBodyPoint>&, int);

	//算法函数
private: 
	//寻找原始nurbs曲线上曲率半径小于给定offset值的点,(如果首位点相同,则记录两次,当作距离为0的线段)
	//所谓ptsCurvRadLower得到的个数一定是偶数
	void FindPtsCurvRadLowerThanOffset(const QVector<NurbsBodyPoint>&, const int&, QVector<int>&);
	//计算凸包
	double CalConvexHull(const QVector<NurbsBodyPoint>&, const int&);
};
