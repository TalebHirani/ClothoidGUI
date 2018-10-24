#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->plot->addGraph();
    ui->plot->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
    ui->plot->graph(0)->setLineStyle(QCPGraph::lsNone);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::addClothoid(double x0, double y0, double theta0, double x1, double y1, double theta1)
{
    double k;
    double dk;
    double L;
    uint npts = 50;
    int res = Clothoid::buildClothoid(x0, y0, theta0, x1, y1, theta1, k, dk, L);
    std::vector<double> X(10), Y(10);
    res = Clothoid::pointsOnClothoid(x0, y0, theta0, k, dk, L, npts, X, Y);
    for(uint i = 0; i < npts; i++) {
        addPoint(X[i],Y[i]);
    }
    plot();
}

void MainWindow::addPoint(double x, double y)
{
    qv_x.append(x);
    qv_y.append(y);
}

void MainWindow::addLine(double x0, double y0, double x1, double y1)
{
    int npts = 20;
    double dy = (y1-y0)/npts;
    double dx = (x1-x0)/npts;
    for(int i = 0; i <= npts; i++) {
        addPoint(x0,y0);
        x0 = x0 + dx;
        y0 = y0 + dy;
    }
    plot();
}

void MainWindow::plot()
{
    ui->plot->graph(0)->setData(qv_x, qv_y);
    ui->plot->yAxis->rescale(true);
    ui->plot->xAxis->rescale(true);
    ui->plot->replot();
    ui->plot->update();
}


void MainWindow::on_btn_Clothoid_clicked()
{

    addClothoid(ui->sb_x0->value(),ui->sb_y0->value(),ui->sb_theta0->value(),ui->sb_x1->value(),ui->sb_y1->value(),ui->sb_theta1->value());
}

void MainWindow::on_btn_Line_clicked()
{
    addLine(ui->sb_x0->value(), ui->sb_y0->value(), ui->sb_x1->value(), ui->sb_y1->value());
    plot();
}

void MainWindow::on_btn_Export_clicked()
{
      QDir::setCurrent("../../../..");
      QString filename = "output.txt";
      QFile file(filename);
      if (file.open(QIODevice::ReadWrite)) {
              QTextStream stream(&file);
              for(int i = 0; i < qv_x.length(); i++) {
                stream << qv_x[i] << "               " << qv_y[i] << endl;
              }

          }
}
