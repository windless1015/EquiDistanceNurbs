#include "EquiDistanceNurbs.h"
#include <qevent.h>
#include <qpainter.h>
#include <qtextstream.h>
#include <qdebug.h>
#include <QIntValidator>

EquiDistanceNurbs::EquiDistanceNurbs(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	isShowCtrlPoints = true;
	isShowBodyPtsTangent = false;
	isShowCtrlPtsConnectedLine = false;
	isShowCurvatureRadius = false;
	isShowBodyPtsNorm = false;

	//lineedit只能输入数字
	ui.lineEdit_idx->setValidator(new QIntValidator(0, 10000, this));

	connect(ui.btnDisplayTan, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayTangent);
	connect(ui.btnDisplayNorm, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayNorm);
	connect(ui.btnOffsetCurve, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnEquidistance);
	connect(ui.btnOffsetCurve, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnEquidistance);
	connect(ui.btnDisplayCurvatureRad, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayCurvature);
	connect(ui.lineEdit_idx, &QLineEdit::returnPressed, this, &EquiDistanceNurbs::OnDisplayCurvature);
}

void EquiDistanceNurbs::OnBtnDisplayCurvature()
{
	isShowCurvatureRadius = !isShowCurvatureRadius;
	update();
}

void EquiDistanceNurbs::OnDisplayCurvature()
{
	int lineEditNum = ui.lineEdit_idx->text().toInt();
	if (lineEditNum < 0 || lineEditNum >= originNurbsBodyPts.size())
		return;
	update();
}

void EquiDistanceNurbs::OnBtnDisplayTangent()
{
	isShowBodyPtsTangent = !isShowBodyPtsTangent;
	update();
}

void EquiDistanceNurbs::OnBtnDisplayNorm()
{
	isShowBodyPtsNorm = !isShowBodyPtsNorm;
	update();
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
		CalFirstDerivativeAndNorm(originNurbsBodyPts);
		CalSecondDerivativeAndCurvRad(originNurbsBodyPts);
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
	curvePen.setWidth(2);
	painter.setPen(curvePen);

	for (int i = 0; i < originNurbsBodyPts.size() - 1; i++) //遍历body pts是从idx为[0, i-2]
	{
		//绘制每一小段的直线,当body pts足够多的时候就是曲线
		painter.drawLine(originNurbsBodyPts[i].point, originNurbsBodyPts[i + 1].point);
	}

	//显示法向量
	if (isShowBodyPtsNorm)
	{
		QPen normPen(QColor(1, 0, 0));
		normPen.setWidth(1);
		painter.setPen(normPen);
		for (int i = 0; i < originNurbsBodyPts.size(); i++)
		{
			QPointF normPt = 100 * originNurbsBodyPts[i].unitNormVector + originNurbsBodyPts[i].point;
			//绘制每一小段的直线,当body pts足够多的时候就是曲线
			painter.drawLine(originNurbsBodyPts[i].point, normPt);
		}
	}

	//显示切向量
	if (isShowBodyPtsTangent)
	{
		QPen tanPen(QColor(1, 0, 0));
		tanPen.setWidth(1);
		//painter.setPen(tanPen);


		painter.setPen(QColor(Qt::blue));
		painter.setBrush(QBrush(Qt::white));
		for (int i = 0; i < originNurbsBodyPts.size(); i++)
		{
			QPointF unitTangent = originNurbsBodyPts[i].firstDerivative;
			unitTangent /= sqrt(unitTangent.x() * unitTangent.x() + unitTangent.y() * unitTangent.y());
			QPointF tangentPt = 80 * unitTangent + originNurbsBodyPts[i].point;
			//绘制每一小段的直线,当body pts足够多的时候就是曲线
			painter.drawLine(originNurbsBodyPts[i].point, tangentPt);
		}

	}

	//显示曲率半径
	if (isShowCurvatureRadius)
	{
		QPen curvaturePen(QColor(0, 0, 0));
		curvaturePen.setWidth(1);
		painter.setPen(curvaturePen);
		painter.setBrush(QBrush(Qt::white));
		int lineEditNum = ui.lineEdit_idx->text().toInt();
		int r = static_cast<int>(originNurbsBodyPts[lineEditNum].curvatureRadius);
		QPointF circleCenter = r * originNurbsBodyPts[lineEditNum].unitNormVector + originNurbsBodyPts[lineEditNum].point;
		painter.drawEllipse(circleCenter, r, r);
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


void EquiDistanceNurbs::CalFirstDerivativeAndNorm(QVector<NurbsBodyPoint>& bodyPtsVector)
{
	int vectorSize = bodyPtsVector.size();
	if (vectorSize <= 0)
		return;
	//从1开始到 vectorSize - 2结束
	for (int i = 1; i < vectorSize - 1; i++)
	{
		double firstDerivative_x = (bodyPtsVector[i + 1].point.x() - bodyPtsVector[i - 1].point.x()) / (2 * step);
		double firstDerivative_y = (bodyPtsVector[i + 1].point.y() - bodyPtsVector[i - 1].point.y()) / (2 * step);
		//计算x和y方向的一阶导数
		bodyPtsVector[i].firstDerivative.rx() = firstDerivative_x;
		bodyPtsVector[i].firstDerivative.ry() = firstDerivative_y;
		//计算法向量, 参考论文1.2章节, 曲线法向
		double theta = sqrt(firstDerivative_x* firstDerivative_x + firstDerivative_y * firstDerivative_y);
		//单位化
		bodyPtsVector[i].unitNormVector.rx() = -firstDerivative_y / theta;
		bodyPtsVector[i].unitNormVector.ry() = firstDerivative_x / theta;
	}

	////再处理曲线的两个端点, 点的坐标不用变化, 其他的信息均可以以相邻点信息替换
	bodyPtsVector[0].firstDerivative = bodyPtsVector[1].firstDerivative;
	bodyPtsVector[0].unitNormVector = bodyPtsVector[1].unitNormVector;
	bodyPtsVector[vectorSize - 1].firstDerivative = bodyPtsVector[vectorSize - 2].firstDerivative;
	bodyPtsVector[vectorSize - 1].unitNormVector = bodyPtsVector[vectorSize - 2].unitNormVector;
}

void EquiDistanceNurbs::CalSecondDerivativeAndCurvRad(QVector<NurbsBodyPoint>&bodyPtsVector)
{
	//计算二阶导数是一样的
	int vectorSize = bodyPtsVector.size();
	if (vectorSize <= 0)
		return;
	for (int i = 1; i < vectorSize - 1; i++)
	{
		double secondDerivative_x = (bodyPtsVector[i + 1].firstDerivative.x() - bodyPtsVector[i - 1].firstDerivative.x()) / (2 * step);
		double secondDerivative_y = (bodyPtsVector[i + 1].firstDerivative.y() - bodyPtsVector[i - 1].firstDerivative.y()) / (2 * step);

		//计算x和y方向的二阶导数
		bodyPtsVector[i].secondDerivative.rx() = secondDerivative_x;
		bodyPtsVector[i].secondDerivative.ry() = secondDerivative_y;
		//计算曲率半径
		double temp1 = bodyPtsVector[i].firstDerivative.x() * bodyPtsVector[i].firstDerivative.x();
		double temp2 = bodyPtsVector[i].firstDerivative.y() * bodyPtsVector[i].firstDerivative.y();
		double num = pow(temp1 + temp2, 1.5);//分子
		double temp3 = bodyPtsVector[i].firstDerivative.x() * bodyPtsVector[i].secondDerivative.y();
		double temp4 = bodyPtsVector[i].secondDerivative.x() * bodyPtsVector[i].firstDerivative.y();
		double den = qAbs(temp3 - temp4); //论文里面的曲率半径的公式给的是错的
		double curvatureRadius;
		curvatureRadius = num / den;
		if (curvatureRadius > DMAX)
			curvatureRadius = DMAX;
		bodyPtsVector[i].curvatureRadius = curvatureRadius;
	}
	bodyPtsVector[0].secondDerivative = bodyPtsVector[1].secondDerivative;
	bodyPtsVector[0].curvatureRadius = bodyPtsVector[1].curvatureRadius;
	bodyPtsVector[vectorSize - 1].secondDerivative = bodyPtsVector[vectorSize - 2].secondDerivative;
	bodyPtsVector[vectorSize - 1].curvatureRadius = bodyPtsVector[vectorSize - 2].curvatureRadius;
}