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
	//���� ����
	const int DEGREE = 3;
	//��֮��Ĳ�������
	const double step = 0.01;
	//ԭʼnurbs���ߵĿ��Ƶ�
	QVector<NurbsCtrlPoint>			nurbsCtrlPoints;
	//ԭʼnurbs�����ϵĵ�
	QVector<NurbsBodyPoint>		originNurbsBodyPts;
	QVector<NurbsBodyPoint> locationVec;

private:
	bool isShowCtrlPoints; //�Ƿ���ʾ���Ƶ�
	bool isShowCtrlPtsConnectedLine; //�Ƿ���ʾ���Ƶ�֮���ֱ����
	bool isShowBodyPtsTangent; //�Ƿ���ʾ���ߵ��������
	bool isShowBodyPtsNorm; //�Ƿ���ʾ���ߵ�ķ�����
	bool isShowCurvatureRadius; //�Ƿ���ʾĳһ��������ʰ뾶
	bool isShowRawEquiDistance; //�Ƿ���ʾԭʼ�Ⱦ���

private:
	void GenerateNURBSCurve();
	double N3Spline(int i, double u);
	//����һ�׵����ͷ�����
	void CalFirstDerivativeAndNorm(QVector<NurbsBodyPoint>&);
	//������׵��������ʰ뾶
	void CalSecondDerivativeAndCurvRad(QVector<NurbsBodyPoint>&);
	//����ƫ������
	void CalOffsetCurve(const QVector<NurbsBodyPoint>&, int);

	//�㷨����
private: 
	//Ѱ��ԭʼnurbs���������ʰ뾶С�ڸ���offsetֵ�ĵ�,(�����λ����ͬ,���¼����,��������Ϊ0���߶�)
	//��νptsCurvRadLower�õ��ĸ���һ����ż��
	void FindPtsCurvRadLowerThanOffset(const QVector<NurbsBodyPoint>&, const int&, QVector<int>&);
	//����͹��
	double CalConvexHull(const QVector<NurbsBodyPoint>&, const int&);
};
