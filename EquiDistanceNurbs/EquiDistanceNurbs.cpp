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
	isShowRawEquiDistance = false;
	//lineeditֻ����������
	ui.lineEdit_idx->setValidator(new QIntValidator(0, 10000, this));

	connect(ui.btnDisplayTan, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayTangent);
	connect(ui.btnDisplayNorm, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayNorm);
	connect(ui.btnShowOriginEquiLines, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnRawEquidistance);
	connect(ui.btnDisplayCurvatureRad, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayCurvature);
	connect(ui.lineEdit_idx, &QLineEdit::returnPressed, this, &EquiDistanceNurbs::OnDisplayCurvature);

	NurbsCtrlPoint p1;
	p1.point = QPointF(243, 676);
	p1.weight = 1.0;
	nurbsCtrlPoints.push_back(p1);
	p1.point = QPointF(244, 304);
	nurbsCtrlPoints.push_back(p1);
	p1.point = QPointF(356, 672);
	nurbsCtrlPoints.push_back(p1);
	p1.point = QPointF(446, 290);
	nurbsCtrlPoints.push_back(p1);
	p1.point = QPointF(450, 695);
	nurbsCtrlPoints.push_back(p1);
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

void EquiDistanceNurbs::OnBtnRawEquidistance()
{

	isShowRawEquiDistance = !isShowRawEquiDistance;
	CalOffsetCurve(originNurbsBodyPts, 40);

	update();
}

void EquiDistanceNurbs::mousePressEvent(QMouseEvent* event)
{
	if (event->button() == Qt::LeftButton)
	{
		NurbsCtrlPoint oneCtrlPt;
		oneCtrlPt.point = event->pos();;
		oneCtrlPt.weight = 1.0;
		//nurbsCtrlPoints.push_back(oneCtrlPt);
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
		QBrush brush(Qt::SolidPattern);//��ˢ
		brush.setColor(Qt::blue);
		painter.setBrush(brush);//���û�ˢ
		for (int i = 0; i < nurbsCtrlPoints.size(); i++)
		{
			//���Ƶ�ʹ��drawEllipse ,��ctrlPtsPen.setWidth����Ϊ5(�����) �����ʾ������
			painter.drawEllipse(static_cast<int>(nurbsCtrlPoints[i].point.x()) ,
				static_cast<int>(nurbsCtrlPoints[i].point.y()) ,2, 2);
			//painter.drawPoint(nurbsCtrlPoints[i].point.x(), nurbsCtrlPoints[i].point.y());
		}
	}

	//��ʾ���Ƶ�֮��ֱ����
	if (isShowCtrlPtsConnectedLine)
	{
		QPen ctrlPen(QColor(0, 255, 0));
		ctrlPen.setWidth(2);
		ctrlPen.setStyle(Qt::DotLine); //����ֱ����ʽΪ������
		painter.setPen(ctrlPen);
		for (int i = 0; i < nurbsCtrlPoints.size() - 1; i++)
		{
			painter.drawLine(nurbsCtrlPoints[i].point, nurbsCtrlPoints[i + 1].point);
		}
	}

	//����NURBS����
	QPen curvePen(QColor(0, 0, 0));
	curvePen.setWidth(2);
	painter.setPen(curvePen);

	for (int i = 0; i < originNurbsBodyPts.size() - 1; i++) //����body pts�Ǵ�idxΪ[0, i-2]
	{
		//����ÿһС�ε�ֱ��,��body pts�㹻���ʱ���������
		painter.drawLine(originNurbsBodyPts[i].point, originNurbsBodyPts[i + 1].point);
	}

	//��ʾ������
	if (isShowBodyPtsNorm)
	{
		QPen normPen(QColor(1, 0, 0));
		normPen.setWidth(1);
		painter.setPen(normPen);
		for (int i = 0; i < originNurbsBodyPts.size(); i++)
		{
			QPointF normPt = 100 * originNurbsBodyPts[i].unitNormVector + originNurbsBodyPts[i].point;
			//����ÿһС�ε�ֱ��,��body pts�㹻���ʱ���������
			painter.drawLine(originNurbsBodyPts[i].point, normPt);
		}
	}

	//��ʾ������
	if (isShowBodyPtsTangent)
	{
		QPen tanPen(QColor(1, 0, 0));
		tanPen.setWidth(1);
		painter.setPen(QColor(Qt::blue));
		painter.setBrush(QBrush(Qt::white));
		QPointF unitTangent;
		for (int i = 0; i < originNurbsBodyPts.size(); i++)
		{
			unitTangent = originNurbsBodyPts[i].firstDerivative;
			unitTangent /= sqrt(unitTangent.x() * unitTangent.x() + unitTangent.y() * unitTangent.y());
			QPointF tangentPt = 80 * unitTangent + originNurbsBodyPts[i].point;
			//����ÿһС�ε�ֱ��,��body pts�㹻���ʱ���������
			painter.drawLine(originNurbsBodyPts[i].point, tangentPt);
		}

	}

	//��ʾ���ʰ뾶
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

	//����ԭʼ��nurb������ʾraw�Ⱦ���(���Խ�,δ�����ι���)
	if (isShowRawEquiDistance)
	{
		QPen rawEquiPen(QColor(1, 0, 0));
		rawEquiPen.setWidth(1);
		painter.setPen(QColor(Qt::blue));
		QPointF rawEquiNurbPt_prev;
		QPointF rawEquiNurbPt_follw;
		int vecSize = originNurbsBodyPts.size();
		for (int i = 0; i < vecSize - 2; i++)
		{
			rawEquiNurbPt_prev = originNurbsBodyPts[i].unitNormVector * 40 + originNurbsBodyPts[i].point;
			rawEquiNurbPt_follw = originNurbsBodyPts[i+1].unitNormVector * 40 + originNurbsBodyPts[i+1].point;
			//����ÿһС�ε�ֱ��,��body pts�㹻���ʱ���������
			painter.drawLine(rawEquiNurbPt_prev, rawEquiNurbPt_follw);
		}
		//�ٻ������һ��
		rawEquiNurbPt_prev = rawEquiNurbPt_follw;
		rawEquiNurbPt_follw = originNurbsBodyPts[vecSize - 1].unitNormVector * 40 + originNurbsBodyPts[vecSize - 1].point;
		painter.drawLine(rawEquiNurbPt_prev, rawEquiNurbPt_follw);
	
	}

	//��ʱ��ʾ
	{
		QPen ctrlPtsPen(QColor(0, 0, 255));
		ctrlPtsPen.setWidth(5);
		painter.setPen(ctrlPtsPen);
		QBrush brush(Qt::SolidPattern);//��ˢ
		brush.setColor(Qt::blue);
		painter.setBrush(brush);//���û�ˢ
		for (int i = 0; i < locationVec.size(); i++)
		{
			QPointF tmpPt = locationVec[i].point;
			tmpPt = tmpPt + locationVec[i].unitNormVector * 40;
			painter.drawPoint(tmpPt.x(), tmpPt.y());
		}
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
	//��1��ʼ�� vectorSize - 2����
	for (int i = 1; i < vectorSize - 1; i++)
	{
		double firstDerivative_x = (bodyPtsVector[i + 1].point.x() - bodyPtsVector[i - 1].point.x()) / (2 * step);
		double firstDerivative_y = (bodyPtsVector[i + 1].point.y() - bodyPtsVector[i - 1].point.y()) / (2 * step);
		//����x��y�����һ�׵���
		bodyPtsVector[i].firstDerivative.rx() = firstDerivative_x;
		bodyPtsVector[i].firstDerivative.ry() = firstDerivative_y;
		//���㷨����, �ο�����1.2�½�, ���߷���
		double theta = sqrt(firstDerivative_x* firstDerivative_x + firstDerivative_y * firstDerivative_y);
		//��λ��
		bodyPtsVector[i].unitNormVector.rx() = -firstDerivative_y / theta;
		bodyPtsVector[i].unitNormVector.ry() = firstDerivative_x / theta;
	}

	////�ٴ������ߵ������˵�, ������겻�ñ仯, ��������Ϣ�����������ڵ���Ϣ�滻
	bodyPtsVector[0].firstDerivative = bodyPtsVector[1].firstDerivative;
	bodyPtsVector[0].unitNormVector = bodyPtsVector[1].unitNormVector;
	bodyPtsVector[vectorSize - 1].firstDerivative = bodyPtsVector[vectorSize - 2].firstDerivative;
	bodyPtsVector[vectorSize - 1].unitNormVector = bodyPtsVector[vectorSize - 2].unitNormVector;
}

void EquiDistanceNurbs::CalSecondDerivativeAndCurvRad(QVector<NurbsBodyPoint>&bodyPtsVector)
{
	//������׵�����һ����
	int vectorSize = bodyPtsVector.size();
	if (vectorSize <= 0)
		return;
	for (int i = 1; i < vectorSize - 1; i++)
	{
		double secondDerivative_x = (bodyPtsVector[i + 1].firstDerivative.x() - bodyPtsVector[i - 1].firstDerivative.x()) / (2 * step);
		double secondDerivative_y = (bodyPtsVector[i + 1].firstDerivative.y() - bodyPtsVector[i - 1].firstDerivative.y()) / (2 * step);

		//����x��y����Ķ��׵���
		bodyPtsVector[i].secondDerivative.rx() = secondDerivative_x;
		bodyPtsVector[i].secondDerivative.ry() = secondDerivative_y;
		//�������ʰ뾶
		double temp1 = bodyPtsVector[i].firstDerivative.x() * bodyPtsVector[i].firstDerivative.x();
		double temp2 = bodyPtsVector[i].firstDerivative.y() * bodyPtsVector[i].firstDerivative.y();
		double num = pow(temp1 + temp2, 1.5);//����
		double temp3 = bodyPtsVector[i].firstDerivative.x() * bodyPtsVector[i].secondDerivative.y();
		double temp4 = bodyPtsVector[i].secondDerivative.x() * bodyPtsVector[i].firstDerivative.y();
		double den = qAbs(temp3 - temp4); //������������ʰ뾶�Ĺ�ʽ�����Ǵ��
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

///////////////////////////////////�㷨����///////////////////////////
void EquiDistanceNurbs::FindPtsCurvRadLowerThanOffset(const QVector<NurbsBodyPoint>&nurbsBodyPts,
	const int& offset, QVector<int>& ptsCurvRadLower)
{
	int nurbsBodyPtSize = nurbsBodyPts.size();
	if (nurbsBodyPtSize <= 0)
		return;
	int j = 0;
	double curvRad;
	while (j < nurbsBodyPtSize)
	{
		while (nurbsBodyPts[j].curvatureRadius > qAbs(offset) && j <= nurbsBodyPtSize - 2) //������������2����
		{
			j++;
		}
		if (j == nurbsBodyPtSize - 1) //���������һ���������
		{
			break;
		}
		ptsCurvRadLower.push_back(j);//���
		while (nurbsBodyPts[j].curvatureRadius <= qAbs(offset) && j < nurbsBodyPtSize)
		{
			j++;
		}
		ptsCurvRadLower.push_back(j - 1); //����һ������, �յ�
	}

}

void EquiDistanceNurbs::CalOffsetCurve(const QVector<NurbsBodyPoint>&bodyPtsVector, int offset)
{
	//1.��ȡnurbs������С�����ʰ뾶����ĩ�˵�
    QVector<int> startEndPtLowerCurvRad;       
	FindPtsCurvRadLowerThanOffset(bodyPtsVector, offset, startEndPtLowerCurvRad);

	//2.��startEndPtLowerCurvRad��ŵ����ʰ뾶С��offset�ĵ�, �޳���ͬ�����purifyEndPts
	QVector<int> purifyEndPts;
	if (startEndPtLowerCurvRad.size() > 0)
	{
		int comparedIdx = startEndPtLowerCurvRad.at(0); //��ȡ��һ�����Ƚϵ�idx
		purifyEndPts.push_back(comparedIdx); //��һ�����ݿ϶�����,ѹ��
		for (int i = 1; i < startEndPtLowerCurvRad.size(); i++)
		{
			int curIdx = startEndPtLowerCurvRad[i];
			if (comparedIdx != curIdx)
			{
				purifyEndPts.push_back(curIdx);
				comparedIdx = curIdx;
			}
		}
	}
	//3. ��purifyEndPts�в������Խ��ĵ�ɾ��,ʣ�µ����Խ��ĵ�
	double con;
	QVector<int> needCutEndPts;  //��Ҫ�ü�����ĩ�˵�
	for (int i = 0; i < purifyEndPts.size(); i++)
	{
		con = CalConvexHull(bodyPtsVector, purifyEndPts.at(i));
		if (offset > 0 && con > 0)
		{
			needCutEndPts.push_back(purifyEndPts.at(i));
		}
		else if (offset < 0 && con < 0)
		{
			needCutEndPts.push_back(purifyEndPts.at(i));
		}
	}


	//��ʾ��ƫ�Ƶ�
	for (int i = 0; i < needCutEndPts.size(); i++)
	{
		locationVec.push_back(bodyPtsVector.at(needCutEndPts.at(i)));
	}

}


double EquiDistanceNurbs::CalConvexHull(const QVector<NurbsBodyPoint>&bodyPtsVector, const int& index)
{
	if (index < 1)
		return 0;
	QPointF p0, p1, p2;
	double result = 0;
	p0 = bodyPtsVector.at(index - 1).point;
	p1 = bodyPtsVector.at(index).point;
	p2 = bodyPtsVector.at(index + 1).point;

	//���ʾ�������������� ���ж��Ƿ���͹�����ǰ���
	double denom1, denom2;
	denom1 = p2.x() - p1.x();
	denom2 = p1.y() - p0.y();
	if ((fabs(denom1) > 1e-5) && (fabs(denom2) > 1e-5))
		result = (p1.x() - p0.x())*(p2.y() - p1.y()) - denom1 * denom2;  //���� ʽ6������ʽ,�ж�������
	return result;
}
