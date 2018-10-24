#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFile>
#include <QTextStream>
#include "clothoid.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void addClothoid(double x0, double y0, double theta0, double x1, double y1, double theta1);
    void addPoint(double x, double y);
    void addLine(double x0, double y0, double x1, double y1);
    void plot();

private slots:

    void on_btn_Clothoid_clicked();

    void on_btn_Line_clicked();

    void on_btn_Export_clicked();

private:
    Ui::MainWindow *ui;

    QVector<double> qv_x, qv_y;
};

#endif // MAINWINDOW_H
