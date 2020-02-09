#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_EquiDistanceNurbs.h"
#include <qvector.h>

//������Сֵ,���ֵ
const double DMAX = 1.e12;
const double DMIN = 1.e-12;

typedef struct WeightPoint
{
	QPointF point;
	double weight;
} NurbsCtrlPoint;

typedef struct BodyPoint
{
	QPointF point; //�������
	QPointF firstDerivative; //һ�׵���, x��ʾx�����һ�׵���, y��ʾy�����һ�׵���
	QPointF secondDerivative; //���׵���,x��ʾx����Ķ��׵���, y��ʾy����Ķ��׵���
	QPointF unitNormVector; //��λ������
	double curvatureRadius; //���ʰ뾶
} NurbsBodyPoint;


class EquiDistanceNurbs : public QMainWindow
{
	Q_OBJECT

public:
	EquiDistanceNurbs(QWidget *parent = Q_NULLPTR);



public slots:
	void OnBtnEquidistance();
	void OnBtnDisplayTangent();
	void OnBtnDisplayNorm();
	void OnBtnDisplayCurvature();
	void OnDisplayCurvature();


protected:
	void mousePressEvent(QMouseEvent* event);
	void paintEvent(QPaintEvent* event);

private:
	Ui::EquiDistanceNurbsClass ui;
	//���� ����
	const int DEGREE = 3;
	//��֮��Ĳ�������
	const double step = 0.1;

	//ԭʼnurbs���ߵĿ��Ƶ�
	QVector<NurbsCtrlPoint>			nurbsCtrlPoints;
	//ԭʼnurbs�����ϵĵ�
	QVector<NurbsBodyPoint>		originNurbsBodyPts;

private:
	bool isShowCtrlPoints; //�Ƿ���ʾ���Ƶ�
	bool isShowCtrlPtsConnectedLine; //�Ƿ���ʾ���Ƶ�֮���ֱ����
	bool isShowBodyPtsTangent; //�Ƿ���ʾ���ߵ��������
	bool isShowBodyPtsNorm; //�Ƿ���ʾ���ߵ�ķ�����
	bool isShowCurvatureRadius; //�Ƿ���ʾĳһ��������ʰ뾶

private:
	void GenerateNURBSCurve();
	double N3Spline(int i, double u);
	//����һ�׵����ͷ�����
	void CalFirstDerivativeAndNorm(QVector<NurbsBodyPoint>&);
	//������׵��������ʰ뾶
	void CalSecondDerivativeAndCurvRad(QVector<NurbsBodyPoint>&);
};
