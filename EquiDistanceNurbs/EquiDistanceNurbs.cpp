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
	//lineedit只能输入数字
	ui.lineEdit_idx->setValidator(new QIntValidator(0, 10000, this));

	connect(ui.btnDisplayTan, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayTangent);
	connect(ui.btnDisplayNorm, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayNorm);
	connect(ui.btnShowOriginEquiLines, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnRawEquidistance);
	connect(ui.btnDisplayCurvatureRad, &QPushButton::clicked, this, &EquiDistanceNurbs::OnBtnDisplayCurvature);
	connect(ui.lineEdit_idx, &QLineEdit::returnPressed, this, &EquiDistanceNurbs::OnDisplayCurvature);

	/*NurbsCtrlPoint p1;
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
	nurbsCtrlPoints.push_back(p1);*/
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
	CalOffsetCurve(originNurbsBodyPts, 40);
	isShowRawEquiDistance = !isShowRawEquiDistance;
	update();
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
		painter.setPen(QColor(Qt::blue));
		painter.setBrush(QBrush(Qt::white));
		QPointF unitTangent;
		for (int i = 0; i < originNurbsBodyPts.size(); i++)
		{
			unitTangent = originNurbsBodyPts[i].firstDerivative;
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

	//根据原始的nurb曲线显示raw等距线(有自交,未经修饰过的)
	if (isShowRawEquiDistance && rawOffsetNurbsBodyPts.size() > 0)
	{
		QPen rawEquiPen(QColor(1, 0, 0));
		rawEquiPen.setWidth(1);
		painter.setPen(QColor(Qt::blue));
		QPointF rawEquiNurbPt_prev;
		QPointF rawEquiNurbPt_follw;
		int vecSize = locationVec.size(); //偏移曲线的数量
		/*for (int i = 0; i < vecSize - 2; i++)
		{
			rawEquiNurbPt_prev = rawOffsetNurbsBodyPts.at(i).point;
			rawEquiNurbPt_follw = rawOffsetNurbsBodyPts.at(i+1).point;
			painter.drawLine(rawEquiNurbPt_prev, rawEquiNurbPt_follw);
		}*/
		for (int i = 0; i < vecSize - 1; i++)
		{
			rawEquiNurbPt_prev = locationVec.at(i).point;
			rawEquiNurbPt_follw = locationVec.at(i+1).point;
			painter.drawLine(rawEquiNurbPt_prev, rawEquiNurbPt_follw);
		}

	}

	//临时显示
	{
		return;
		QPen ctrlPtsPen(QColor(0, 0, 255));
		ctrlPtsPen.setWidth(5);
		painter.setPen(ctrlPtsPen);
		QBrush brush(Qt::SolidPattern);//画刷
		brush.setColor(Qt::blue);
		painter.setBrush(brush);//设置画刷
		for (int i = 0; i < locationVec.size(); i++)
		{
			//QPointF tmpPt = locationVec[i].point;
			//tmpPt = tmpPt + locationVec[i].unitNormVector * 40;
			//painter.drawPoint(tmpPt.x(), tmpPt.y());
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

///////////////////////////////////算法函数///////////////////////////
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
		while (nurbsBodyPts[j].curvatureRadius > qAbs(offset) && j <= nurbsBodyPtSize - 2) //遍历到导数第2个点
		{
			j++;
		}
		if (j == nurbsBodyPtSize - 1) //遍历到最后一个点就跳出
		{
			break;
		}
		ptsCurvRadLower.push_back(j);//起点
		while (nurbsBodyPts[j].curvatureRadius <= qAbs(offset) && j < nurbsBodyPtSize)
		{
			j++;
		}
		ptsCurvRadLower.push_back(j - 1); //回退一个索引, 终点
	}

}

void EquiDistanceNurbs::CalOffsetCurve(const QVector<NurbsBodyPoint>&bodyPtsVector, int offset)
{

	//先计算原始的nurbs偏移曲线
	rawOffsetNurbsBodyPts.clear();
	NurbsBodyPoint offsetNurbsBodyPt;
	for (int i=0; i < bodyPtsVector.size(); i++)
	{
		offsetNurbsBodyPt.point = originNurbsBodyPts[i].unitNormVector * 40 + originNurbsBodyPts[i].point;
		offsetNurbsBodyPt.unitNormVector = originNurbsBodyPts[i].unitNormVector * -1;//反向
		rawOffsetNurbsBodyPts.push_back(offsetNurbsBodyPt);
	}

	//1.获取nurbs曲线上小于曲率半径的首末端点
    QVector<int> startEndPtLowerCurvRad;       
	FindPtsCurvRadLowerThanOffset(bodyPtsVector, offset, startEndPtLowerCurvRad);

	//2.将startEndPtLowerCurvRad存放的曲率半径小于offset的点, 剔除相同点放入purifyEndPts
	QVector<int> purifyEndPts;
	if (startEndPtLowerCurvRad.size() > 0)
	{
		int comparedIdx = startEndPtLowerCurvRad.at(0); //获取第一个待比较的idx
		purifyEndPts.push_back(comparedIdx); //第一个数据肯定存在,压入
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
	//3. 将purifyEndPts中不产生自交的点删除,剩下的是自交三角形的两个腰角的顶点
	//原理: 根据凸包性判断, 有一些曲率半径比offset小,但是因为是某种凹凸性(可能是凸,也可能是凹),这些点不会产生自交
	double con;
	QVector<int> needCutEndPts;  //需要裁剪的首末端点
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

	//4.寻找自交点
	QVector<int> selfCrossPts; //自交点的idx, 每两个为一对, 分别为大小索引
	CalSelfCrossPoint(rawOffsetNurbsBodyPts, needCutEndPts, selfCrossPts);
	QMap<int, QVector<QPointF>> connectCurves;
	ConstructConnectCurve(rawOffsetNurbsBodyPts, selfCrossPts, connectCurves);
	//5.创建新的经过光顺的偏置曲线
	QVector<NurbsBodyPoint> finalOffsetCurve;
	ConstructCurvesSnippets(connectCurves, selfCrossPts, rawOffsetNurbsBodyPts, finalOffsetCurve);

	//显示其偏移点
	for (int i = 0; i < finalOffsetCurve.size(); i++)
	{
		locationVec.push_back(finalOffsetCurve.at(i));
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

	//本质就是利用向量叉乘 来判断是否是凸包还是凹的
	double denom1, denom2;
	denom1 = p2.x() - p1.x();
	denom2 = p1.y() - p0.y();
	if ((fabs(denom1) > 1e-5) && (fabs(denom2) > 1e-5))
		result = (p1.x() - p0.x())*(p2.y() - p1.y()) - denom1 * denom2;  //论文 式6的行列式,判断正负用
	return result;
}

/*
1.先说作者的方法,作者的方法是找到自交倒三角形的腰点,然后回溯. 使用了三重循环.计算两个点之间距离,设定阈值判断自交点
*/
void EquiDistanceNurbs::CalSelfCrossPoint(const QVector<NurbsBodyPoint>&offsetCurvBodyPts, QVector<int>& needCutEndPts,
	QVector<int>&selfCrossPtVector)
{
	//向只包含腰点坐标数组插入两个idx, 一个是0, 一个是rawOffsetNurbsBodyPts的数组长度
	//因为这样可以知道前后的边界处
	needCutEndPts.push_front(0);
	needCutEndPts.push_back(offsetCurvBodyPts.size());

	//staEndPtIdx表示从这一组腰点中开始遍历, 并且从1开始遍历,因为第一个数据为0(起始边界)
	//并且,遍历的时候需要跨越一组(两个)点, 后续可以优化成一个数据结构
	//遍历到needCutEndPts.size() - 1 最后一个数据是结束边界, 不需要遍历
	double deltaX, deltaY, distance;
	for (int k = 1; k < needCutEndPts.size() - 1; k += 2) 
	{
		bool isFindCrossPt = false; //是否找到自交点标记位
		//最多减少到前一组数据的下边界(因为不会跨越上一个idx较大的腰点) 和 1/step(100次)
		for (int j = 0, lowerPtIdx = needCutEndPts[k];
			needCutEndPts[k - 1] < lowerPtIdx && j < 1 / step; j++, lowerPtIdx--)
		{
			for (int i = 0, greaterPtIdx = needCutEndPts[k + 1];
				greaterPtIdx < needCutEndPts[k + 2] && i < 1 / step; i++, greaterPtIdx++) //同样不会超越下一组数据的上边界
			{
				deltaX = offsetCurvBodyPts[greaterPtIdx].point.x() - offsetCurvBodyPts[lowerPtIdx].point.x();
				deltaY = offsetCurvBodyPts[greaterPtIdx].point.y() - offsetCurvBodyPts[lowerPtIdx].point.y();
				distance = sqrt(deltaX * deltaX + deltaY * deltaY);
				if (distance < 0.7)//距离阈值
				{
					selfCrossPtVector.push_back(lowerPtIdx); //把自交点的两个坐标记录,一个是idx较小的, 一个是较大的,两个点相距非常近,可以认为是一个点
					selfCrossPtVector.push_back(greaterPtIdx);
					isFindCrossPt = true; //找到自交点
					break; //这个是跳出最内层循环
				}
				if (isFindCrossPt) //如果找到一个自交点,那么就要跳出这一次寻找, 寻找下一组点, for循环回到最外层的遍历
				{
					break;
				}
			}
		}
	}
	
}

//rawOffsetPts 表示原始的偏置曲线, 
//crossPts表示相交点的id 值 数组, 一般为偶数
void EquiDistanceNurbs::ConstructConnectCurve(const QVector<NurbsBodyPoint>& rawOffsetPts, 
	QVector<int>&selfCrossPts, QMap<int, QVector<QPointF>>&connectCurves)
{
	if (selfCrossPts.size() <= 0)
		return;
	int deltaBackup = 20;  //后退的步数
	int deltaForward = 20; //向前的步数
	selfCrossPts.push_front(0);
	selfCrossPts.push_back(rawOffsetPts.size());
	//本质就是用一段贝塞尔曲线来做拟合, 这一段曲线能和原曲线连续,1阶,2阶可导. 本质就是先找6个控制点.
	//分别为p0, p1, p2, p3, p4, p5.   其中 p0, p5是连接光顺曲线和原来nurbs曲线的连接点, 其他几个点的求法
	//在论文中均有描述,
	double selfDefRatio = 0.06;       //自定义系数
	int mapCount = 0;
	QPointF p0, p1, p2, p3, p4, p5;  //光顺连接曲线的6个控制点
	//下面带_v的表示向量,不是一个点
	QPointF k_p1_v; //p1点的斜率 (p0 和p1是通过y = kx +b 形式连接,)
	QPointF k_p4_v; //p4点的斜率(p4 和p5是通过y = -kx +b 形式连接,)
	QPointF prevConnPtDeriv_1order_v; //连接点前一个(idx较小的)点的1阶导数
	QPointF connPtDeriv_1order_v; // 连接点的一阶导数
	QPointF followConnPtDeriv_1order_v; //连接点后一个(idx较大)点的一阶导数
	QPointF connPtDeriv_2order_v;//中间点的二阶导数
	//第一个数据是0, 后面每两个是一对, 最后一个是偏移曲线的长度.  例子: 0    (131, 323)    (266,472)     601
	for (int i=0; i< selfCrossPts.size() - 2; i += 2)
	{
		//创建当前 连接的光顺曲线段 
		QVector<QPointF> curConnectCurvePts;
		//真正的nurb曲线和构造光顺曲线段相互连接点
		int realConnPt1 = selfCrossPts[i + 1] - deltaBackup;  //真正的连接点是,把自交点往idx较小方向偏移了一段距离
		int realConnPt2 = selfCrossPts[i + 2] + deltaForward; //向 idx增大方向延伸一段
		selfCrossPts[i + 1] = realConnPt1;   //修改真正的接触点的数据
		selfCrossPts[i + 2] = realConnPt2;
		if (realConnPt1 < 2)// 超过边界点就用边界点来代替(至少保证前面有2点偏移,为了求导数用)
			realConnPt1 = 2;
		if (realConnPt2 > rawOffsetPts.size() - 3)// 超过边界点就用边界点来代替(至少保证后面有2个点)
			realConnPt2 = rawOffsetPts.size() - 3;
		//此时找到了真正的连接点p0和p5
		p0 = rawOffsetPts[realConnPt1].point;
		p5 = rawOffsetPts[realConnPt2].point;

		//下面寻找p1和p2, 这两个点和p0与p5相临近,先计算相关点的导数
		prevConnPtDeriv_1order_v = (p0 - rawOffsetPts[realConnPt1 - 2].point) / (2*step); //连接点前一个点的一阶导数
		connPtDeriv_1order_v = (rawOffsetPts[realConnPt1 +1].point - rawOffsetPts[realConnPt1 - 1].point) / (2 * step);//连接点一阶导
		followConnPtDeriv_1order_v = (rawOffsetPts[realConnPt1 + 2].point - p0) / (2 * step); //连接点前一个点的一阶导数
		connPtDeriv_2order_v = (followConnPtDeriv_1order_v - prevConnPtDeriv_1order_v) / (2 * step); //中间点二阶导
		//计算p1
		k_p1_v = connPtDeriv_1order_v; //p1点是通过p0点由 y = kx +b找到的
		p1 = k_p1_v * selfDefRatio + p0; //这里selfDefRatio 相当于自变量x, 这里设定的是0.06
		//计算p2, //公式在文章1.4 的(4) 有公式
		//bP2.ry() = derivate2_1y * ratio_b*ratio_b / 20 + 2 * bP1.y() - bP0.y();
		p2 = connPtDeriv_2order_v * selfDefRatio * selfDefRatio / 20 + 2 * p1 - p0;
		//下面求p4和p3
		prevConnPtDeriv_1order_v = (p0 - rawOffsetPts[realConnPt2 - 2].point) / (2 * step); //连接点前一个点的一阶导数
		connPtDeriv_1order_v = (rawOffsetPts[realConnPt2 + 1].point - rawOffsetPts[realConnPt2 - 1].point) / (2 * step);//连接点一阶导
		followConnPtDeriv_1order_v = (rawOffsetPts[realConnPt2 + 2].point - p0) / (2 * step); //连接点前一个点的一阶导数
		connPtDeriv_2order_v = (followConnPtDeriv_1order_v - prevConnPtDeriv_1order_v) / (2 * step); //中间点二阶导
		//计算p4
		k_p4_v = connPtDeriv_1order_v;
		p4 = -k_p4_v * selfDefRatio + p5; //这里是p5 - k *x;
		//计算p3
		p3 = connPtDeriv_2order_v * selfDefRatio * selfDefRatio / 20 + 2 * p4 - p5;
		//构建6阶贝塞尔曲线
		QPointF bezierCurvePt;
		for (double t = 0; t <= 1; t += 0.01)
		{
			bezierCurvePt = pow(1 - t, 5)*p0 + 5 * t*pow(1 - t, 4)*p1 + 10 * pow(t, 2)*pow(1 - t, 3)*p2
				+ 10 * pow(t, 3)*pow(1 - t, 2)*p3 + 5 * pow(t, 4)*(1 - t)*p4 + pow(t, 5)*p5;
			curConnectCurvePts.push_back(bezierCurvePt);
		}
		connectCurves.insert(mapCount, curConnectCurvePts);
		mapCount++;
	}


}

void EquiDistanceNurbs::ConstructCurvesSnippets(const QMap<int, QVector<QPointF>>&connectCurves, 
	const QVector<int>&selfCrossPts,const  QVector<NurbsBodyPoint>&rawOffsetPts,
	QVector<NurbsBodyPoint>& outPutOffsetPts)
{
	if (selfCrossPts.size() <= 0)
		return;
	//拼接所有曲线段
	//selfCrossPts 存储的是端点 和  真正接触点
	NurbsBodyPoint nbPt;
	int i = 0;
	for (; i < selfCrossPts.size() - 2; i += 2) //
	{
		for (int snippetIdx = selfCrossPts.at(i); snippetIdx < selfCrossPts.at(i + 1); snippetIdx++)
		{
			nbPt = rawOffsetPts.at(snippetIdx);
			outPutOffsetPts.push_back(nbPt);
		}
		//然后把构造的曲线找出来拼接
		for (int conCurvIdx = 0; conCurvIdx < connectCurves[i / 2].size(); conCurvIdx++)
		{
			nbPt.point = connectCurves[i / 2].at(conCurvIdx);
			outPutOffsetPts.push_back(nbPt);
		}
	}
	//这时候i已经自增到导数第二个数据的idx
	//最后一段
	for (int count = selfCrossPts.at(i); count < selfCrossPts.at(i+1); count++)
	{
		nbPt = rawOffsetPts.at(count);
		outPutOffsetPts.push_back(nbPt);
	}
}