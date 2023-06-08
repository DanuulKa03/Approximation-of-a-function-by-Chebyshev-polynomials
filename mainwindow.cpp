#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{

    resize(MAXSIZEWIDTH, MAXSIZEHEIGHT);
    createGraphs();
    initializationObjectClass();

}

MainWindow::~MainWindow()
{
}


void MainWindow::initializationObjectClass()
{
    setCentralWidget(new QWidget);

    mainLayout = new QHBoxLayout(centralWidget());
    rightLayout = new QVBoxLayout(this);
    additionalRightLayout = new QGridLayout(this);

    QWidget* widget = new QWidget(this);
    widget->setMaximumWidth(this->width()*0.25);
    widget->setLayout(rightLayout);

    QRegExp regex("[-+]?[0]?[.]?[0-9]+");
    QRegExpValidator *validator = new QRegExpValidator(regex, this);

    degreePolynomial = new QLineEdit(this);
    stepLineEdit = new QLineEdit("0.1", this);
    fromLineEdit = new QLineEdit("-1",this);
    upLineEdit = new QLineEdit("1",this);

    degreePolynomial ->setValidator(validator);
    stepLineEdit ->setValidator(validator);
    fromLineEdit->setValidator(validator);
    upLineEdit ->setValidator(validator);

    startGraphsButton = new QPushButton("Start graphs",this);
//    plusButton = new QPushButton("+", this);
//    minusButton = new QPushButton("-", this);

    QLabel* degree = new QLabel("Degree:", this);
    QLabel* step = new QLabel("Step", this);
    QLabel* valueFrom = new QLabel("From:", this);
    QLabel* valueUp = new QLabel("Up:", this);

    mainLayout->addWidget(customPlot);
    mainLayout->addWidget(widget);
//    mainLayout->addWidget(plusButton, 0, Qt::AlignTop | Qt::AlignLeft);
//    mainLayout->addWidget(minusButton, 0, Qt::AlignTop | Qt::AlignLeft);

    additionalRightLayout->addWidget(valueFrom, 0, 0);
    additionalRightLayout->addWidget(valueUp, 0, 1);
    additionalRightLayout->addWidget(fromLineEdit, 1, 0);
    additionalRightLayout->addWidget(upLineEdit, 1, 1);

    rightLayout->addWidget(degree);
    rightLayout->addWidget(degreePolynomial);
    rightLayout->addWidget(step);
    rightLayout->addWidget(stepLineEdit);
    rightLayout->addLayout(additionalRightLayout);
    rightLayout->addWidget(startGraphsButton);
    rightLayout->setAlignment(Qt::AlignTop);

    QObject::connect(startGraphsButton, &QPushButton::clicked, this, &MainWindow::createPolynomial);
}

void MainWindow::createGraphs()
{
    customPlot = new QCustomPlot(this);
    customPlot->xAxis->setRange(-1,1);
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");
    customPlot->yAxis->setRange(-20, 20);

    connect(customPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->xAxis2, SLOT(setRange(QCPRange)));
    connect(customPlot->yAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->yAxis2, SLOT(setRange(QCPRange)));

    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

}

double MainWindow::computedFunction(double x) //функция, которую мы хотим вычислить
{
    return ( 2*pow(x, 2) * std::cos(x) );
//    return (double)std::cos(x);
}


double MainWindow:: polynomialChebyshev(int n, double x) //вычисление полинома
{
    int a = fromLineEdit->text().toInt(), b = upLineEdit->text().toInt();

    double result = (double)std::cos(n * std::acos( (2.0*x - ( b+a )) / (b-a) ));

    return result;
}

double MainWindow::approximation(double x) //апроксимируемая функция
{
//    int a = fromLineEdit->text().toInt(), b = upLineEdit->text().toInt();

    double sum = 0;
    for (int j = 1; j < degreePolynomial->text().toInt(); j++)
    {
        sum += (double)( 2.0/std::acos(-1) ) * simpsonIntegral(0, std::acos(-1), j ) * polynomialChebyshev(j, x);

//        qDebug() << "         " << "chebyshevCoefficient(j)" << qSetRealNumberPrecision( 10 ) << " = " << chebyshevCoefficient(j) << "\n";
//        qDebug() << "         " << "polynomialChebyshev(j, x)" << qSetRealNumberPrecision( 10 ) << " = " << polynomialChebyshev(j, x) << "\n";
//        qDebug() << "         " << "sum_" << j << qSetRealNumberPrecision( 10 ) << " = " << sum << "\n";
    }


    return ( ( 1.0/std::acos(-1) ) * simpsonIntegral(0, std::acos(-1), 0) + sum );
}

double MainWindow::chebyshevNodes(int k)
{
//    qDebug() << "         " << "(2*k - 1)" << " = " << (2*k - 1) << "\n";
//    qDebug() << "         " << "(2 * degreePolynomial->text().toInt()" << " = "<< qSetRealNumberPrecision( 10 ) << (2 * degreePolynomial->text().toInt()) << "\n";
//    qDebug() << "         " << "std::acos(-1)" << " = "<< qSetRealNumberPrecision( 10 ) << std::acos(-1) << "\n";
//    qDebug() << "x_" << qSetRealNumberPrecision( 10 ) << k << " = " << (double)std::cos( ((2*k - 1)* std::acos(-1)) / (2 * degreePolynomial->text().toInt())) << "\n";

    if (fromLineEdit->text().toInt() != -1 && upLineEdit->text().toInt() != 1)
    {
        int a = fromLineEdit->text().toInt(), b = upLineEdit->text().toInt();

        qDebug() << "x_" << qSetRealNumberPrecision( 10 ) << k << " = " << (double)( (a+b)/2 + ((b-a)*std::cos( ((2*k+1)*std::acos(-1))/ (2*degreePolynomial->text().toInt()))/2 ) ) << "\n";

        return (double)( (a+b)/2.0 + ((b-a)*std::cos( ((2.0*k+1)*std::acos(-1)) / ( 2.0*degreePolynomial->text().toInt()))/2.0 ) );
    }
    else
    {
        double temp = (double)std::cos( ((2.0*k - 1)* std::acos(-1)) / (2.0 * degreePolynomial->text().toInt()));

        return temp;
    }
}

double MainWindow::chebyshevCoefficient(int j)
{
    double sum = 0;
    for (int k = 1; k < degreePolynomial->text().toInt(); k++)
    {
        sum += computedFunction( chebyshevNodes(k) ) * polynomialChebyshev(j, chebyshevNodes(k) );
    }
//    qDebug() << "a_" << j << " = " << ( ((double)(2*sum) /degreePolynomial->text().toInt() ) ) << "\n";
    return ( ((double)(2*sum) /degreePolynomial->text().toInt() ) );
}

double MainWindow::function(double x, int j)
{
    if (j==0) return computedFunction(std::cos(x));

    else return computedFunction( std::cos(x) ) * std::cos( j*x );
}

double MainWindow::simpsonIntegral(double a, double b, int j) {
    const size_t n = (upLineEdit->text().toInt() - fromLineEdit->text().toInt()) / abs(stepLineEdit->text().toDouble()); //вычисляем количество N

    const double width = (b-a)/n;

    double simpson_integral = 0;
    for(size_t step = 0; step < n; step++) {
        const double x1 = a + step*width;
        const double x2 = a + (step+1)*width;

        simpson_integral += (x2-x1)/6.0*(function(x1, j) + 4.0*function(0.5*(x1+x2), j) + function(x2, j));
    }

    return simpson_integral;
}

void MainWindow::createPolynomial()
{
    qDebug() << "Проверка" << "\n";

    checkClearPolynomial(); //проверка на наличие элементов в полиноме

    const size_t N = (upLineEdit->text().toInt() - fromLineEdit->text().toInt()) / abs(stepLineEdit->text().toDouble()); //вычисляем количество N

    for (size_t i = 0; i <= N; i++)
    {
        double value = approximation(fromLineEdit->text().toDouble()+i*stepLineEdit->text().toDouble());
        x.append(fromLineEdit->text().toDouble()+i*stepLineEdit->text().toDouble());
        y.append(value);
    }

    for (size_t i = 0; i <= N; i++)
    {
        double value = computedFunction(fromLineEdit->text().toDouble()+i*stepLineEdit->text().toDouble());
        x2.append(fromLineEdit->text().toDouble()+i*stepLineEdit->text().toDouble());
        y2.append(value);
    }

    customPlot->addGraph();
    customPlot->graph(0)->setPen(QPen(Qt::blue));
    customPlot->addGraph();
    customPlot->graph(1)->setPen(QPen(Qt::red));

    customPlot->graph(0)->setData(x, y);
    customPlot->graph(1)->setData(x2, y2);

    customPlot->graph(0)->rescaleAxes();
    customPlot->graph(1)->rescaleAxes();

    customPlot->yAxis->setRange(-1,1);
    customPlot->xAxis->setRange(-1.25,1.25);

    customPlot->replot();

//    qDebug() << "matrixOdds: " << matrixOdds << "\n";
    qDebug() << "Mx:"<< x << "\n";
    qDebug() << "My:"<< y << "\n";

}

void MainWindow::checkClearPolynomial()
{
    if (!x.isEmpty() && !y.isEmpty())
    {
        x.clear();
        y.clear();
        x2.clear();
        y2.clear();
    }
}
