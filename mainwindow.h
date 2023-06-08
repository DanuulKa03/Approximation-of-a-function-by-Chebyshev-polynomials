#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "qcustomplot.h"

#include <QLineEdit>
#include <QMainWindow>
#include <QPushButton>
#include <cmath>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
public slots:
    void createPolynomial();
private:

    const int MAXSIZEWIDTH = 1000, MAXSIZEHEIGHT = 800;
    int N; //количество точек

    QVector<double> x,y, x2, y2;

    QCustomPlot* customPlot = nullptr;

    QLineEdit* degreePolynomial = nullptr;
    QLineEdit* stepLineEdit = nullptr;
    QLineEdit* fromLineEdit = nullptr;
    QLineEdit* upLineEdit = nullptr;

    QTextEdit* strFunction = nullptr;

    QPushButton* startGraphsButton = nullptr;
    QPushButton* plusButton = nullptr;
    QPushButton* minusButton = nullptr;

    QHBoxLayout* mainLayout = nullptr;
    QVBoxLayout* rightLayout = nullptr;
    QGridLayout* additionalRightLayout = nullptr;

    void createGraphs();
    void initializationObjectClass();
    double polynomialChebyshev(int n, double x);
    void checkClearPolynomial();
    void readingFunctionFromText();
    double computedFunction(double x);
    double approximation(double x);
    double chebyshevNodes(int j);
    double chebyshevCoefficient(int j);
    double simpsonIntegral(double a, double b, int j);
    double function(double x, int j);
};
#endif // MAINWINDOW_H
