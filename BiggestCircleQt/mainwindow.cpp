#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <fstream>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QPointF>
#include <QPolygonF>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    scene = new QGraphicsScene(0, 0, 800, 800, this);
    ui->graphicsView->setScene(scene);
    ui->toolBar->addAction(QStringLiteral("Load from file"), this, SLOT(loadFromFile()));
}

MainWindow::~MainWindow()
{
    delete ui;
    if (scene)
        delete scene;
}

void MainWindow::loadFromFile()
{
    QString filename(QFileDialog::getOpenFileName(this, QStringLiteral("Open File"), QString(), QStringLiteral("Text files(*.txt)")));
    if (filename.isNull())
    {
        qDebug() << "No file";
        return;
    }

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::critical(this, QStringLiteral("Error"), QStringLiteral("Could not open file"));
        return;
    }

    QTextStream in(&file);
    qreal x_c = 0, y_c = 0, r = 0;
    in >> x_c >> y_c >> r;
    int n = 0;
    in >> n;
    QPolygonF poly(n);
    for (auto &p : poly)
        in >> p.rx() >> p.ry();
    file.close();

    scene->clear();
    scene->addPolygon(poly);
    scene->addEllipse(x_c - r, y_c - r, 2 * r, 2 * r);
}
