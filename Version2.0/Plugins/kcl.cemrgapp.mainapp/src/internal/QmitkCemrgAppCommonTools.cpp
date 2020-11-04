/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
/*=========================================================================
 *
 * CemrgApp Common Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbench.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkCoreObjectFactory.h>
#include <mitkProgressBar.h>
#include <QmitkIOUtil.h>

// CemrgAppModule
#include <CemrgCommonUtils.h>
#include "QmitkCemrgAppCommonTools.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QSignalMapper>

const std::string QmitkCemrgAppCommonTools::VIEW_ID = "org.mitk.views.cemrgappcommontools";

void QmitkCemrgAppCommonTools::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::LoadMesh);
    connect(m_Controls.button_2, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertToCarto);
    connect(m_Controls.button_3, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertCarpToVtk);
}

void QmitkCemrgAppCommonTools::SetFocus() {
}

void QmitkCemrgAppCommonTools::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void QmitkCemrgAppCommonTools::LoadMesh() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());
    CemrgCommonUtils::AddToStorage(
                CemrgCommonUtils::LoadVTKMesh(path.toStdString()), "Mesh", this->GetDataStorage());
}

void QmitkCemrgAppCommonTools::ConvertToCarto() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());

    if (path.isEmpty() || !path.endsWith(".vtk")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input File!");
        return;
    }

    //Find thresholds file
    double fileMeanBP = 0.0;
    double fileStdvBP = 0.0;
    bool threshFileExist = false;
    QFileInfo fullPathInfo(path);
    if (fullPathInfo.dir().exists("prodThresholds.txt")) {

        //Threshold file
        QString tPath = fullPathInfo.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "prodThresholds.txt";
        ifstream prodFileRead(tPath.toStdString());
        if (prodFileRead.is_open()) {

            std::string line;
            std::vector<std::string> lines;
            while (getline(prodFileRead, line))
                lines.push_back(line);
            fileMeanBP = QString::fromStdString(lines.at(2)).toDouble();
            fileStdvBP = QString::fromStdString(lines.at(3)).toDouble();
            prodFileRead.close();
            threshFileExist = true;

        }//_if
    }//_if

    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    m_CartoUIThresholding.setupUi(inputs);
    connect(m_CartoUIThresholding.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_CartoUIThresholding.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_CartoUIThresholding.radioButton_1, SIGNAL(toggled(bool)), this, SLOT(ConvertToCartoUITextUpdate()));
    connect(m_CartoUIThresholding.comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(ConvertToCartoUIUpdate()));
    m_CartoUIThresholding.lineEdit_2->setPlaceholderText(QString::number(fileMeanBP));
    m_CartoUIThresholding.lineEdit_3->setPlaceholderText(QString::number(fileStdvBP));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        //Methods
        int methodType = m_CartoUIThresholding.radioButton_1->isChecked() ? 1 : 2;
        bool discreteScheme = m_CartoUIThresholding.comboBox->currentIndex() == 0 ? true : false;

        //Thresholds
        bool ok0;
        std::vector<double> thresholds;
        QRegExp separator("(\\ |\\,|\\;|\\:|\\t)");
        QStringList thresholdsInput = m_CartoUIThresholding.lineEdit_1->text().trimmed().split(separator);
        for (QString item:thresholdsInput)
            if (item.isEmpty())
                thresholdsInput.removeOne(item);
        if (discreteScheme && thresholdsInput.size()==0) {
            QMessageBox::warning(NULL, "Attention", "Reverting to default threshold values!");
            thresholds.push_back(methodType == 1 ? 1.20 : 3.0);
            thresholds.push_back(methodType == 1 ? 1.32 : 4.0);
        } else if (discreteScheme && thresholdsInput.count()>2) {
            QMessageBox::warning(NULL, "Attention", "Parsing thresholds failed!\nReverting to default values.");
            thresholds.push_back(methodType == 1 ? 1.20 : 3.0);
            thresholds.push_back(methodType == 1 ? 1.32 : 4.0);
        } else {
            for (QString item:thresholdsInput) {
                double thresh = item.toDouble(&ok0);
                if (!ok0) {
                    QMessageBox::warning(NULL, "Attention", "Parsing thresholds failed!\nReverting to default values.");
                    thresholds.clear();
                    thresholds.push_back(methodType == 1 ? 1.20 : 3.0);
                    thresholds.push_back(methodType == 1 ? 1.32 : 4.0);
                    break;
                } else
                    thresholds.push_back(thresh);
            }//_for
        }//_if

        //BP values
        bool ok1, ok2;
        double meanBP = m_CartoUIThresholding.lineEdit_2->text().toDouble(&ok1);
        double stdvBP = m_CartoUIThresholding.lineEdit_3->text().toDouble(&ok2);
        if (discreteScheme && (!ok1 || !ok2)) {
            if (threshFileExist) {
                meanBP = fileMeanBP;
                stdvBP = fileStdvBP;
                QMessageBox::information(NULL, "Attention", "Blood pool intensity values from the file was successfully restored.");
            } else {
                QMessageBox::warning(NULL, "Attention", "Parsing blood pool intensity values failed! Try again please.");
                return;
            }//_if
        }//_if

        //Conversion
        CemrgCommonUtils::ConvertToCarto(path.toStdString(), thresholds, meanBP, stdvBP, methodType, discreteScheme);
        QMessageBox::information(NULL, "Attention", "Conversion Completed!");
        inputs->deleteLater();

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if
}

void QmitkCemrgAppCommonTools::ConvertToCartoUIUpdate() {

    if (m_CartoUIThresholding.comboBox->currentIndex() == 1) {
        m_CartoUIThresholding.lineEdit_1->setEnabled(false);
        m_CartoUIThresholding.lineEdit_2->setEnabled(false);
        m_CartoUIThresholding.lineEdit_3->setEnabled(false);
        m_CartoUIThresholding.radioButton_1->setEnabled(false);
        m_CartoUIThresholding.radioButton_2->setEnabled(false);
    } else {
        m_CartoUIThresholding.lineEdit_1->setEnabled(true);
        m_CartoUIThresholding.lineEdit_2->setEnabled(true);
        m_CartoUIThresholding.lineEdit_3->setEnabled(true);
        m_CartoUIThresholding.radioButton_1->setEnabled(true);
        m_CartoUIThresholding.radioButton_2->setEnabled(true);
    }//_if
}

void QmitkCemrgAppCommonTools::ConvertToCartoUITextUpdate() {

    if (m_CartoUIThresholding.radioButton_1->isChecked())
        m_CartoUIThresholding.lineEdit_1->setPlaceholderText("1.2; 1.32");
    else
        m_CartoUIThresholding.lineEdit_1->setPlaceholderText("3; 4");
}

void QmitkCemrgAppCommonTools::ConvertCarpToVtk(){
    QString pathElem = "";
    QString pathPts = "";
    pathElem = QFileDialog::getOpenFileName(NULL, "Open Mesh .elem File");
    if (pathElem.isEmpty() || !pathElem.endsWith(".elem")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.elem) File!");
        return;
    }
    QFileInfo fi(pathElem);
    QString dir = fi.absolutePath();
    QString vtkPath = dir + mitk::IOUtil::GetDirectorySeparator() +  fi.baseName() + ".vtk";

    pathPts = QFileDialog::getOpenFileName(NULL, "Open Mesh .pts File", dir.toStdString().c_str());

    if (pathPts.isEmpty() || !pathPts.endsWith(".pts")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.pts) File!");
        return;
    }

    int regionScalarsReply = QMessageBox::question(NULL, "Question",
            "Include region as (cell) scalar field?", QMessageBox::Yes, QMessageBox::No);

    CemrgCommonUtils::CarpToVtk(pathElem, pathPts, vtkPath, (regionScalarsReply==QMessageBox::Yes));

    int appendScalarFieldReply = QMessageBox::question(NULL, "Question",
            "Append a scalar field from a file?", QMessageBox::Yes, QMessageBox::No);

    if (appendScalarFieldReply==QMessageBox::Yes){
        QString path="";
        QString typeData="";
        int nElem = CemrgCommonUtils::GetTotalFromCarpFile(pathElem);
        int nPts = CemrgCommonUtils::GetTotalFromCarpFile(pathPts);
        int nField;

        while (appendScalarFieldReply==QMessageBox::Yes){
            path = QFileDialog::getOpenFileName(NULL, "Open Scalar field (.dat) file", dir.toStdString().c_str());
            QFileInfo fi2(path);
            std::vector<double> field = CemrgCommonUtils::ReadScalarField(path);

            nField = field.size();
            MITK_INFO << ("FieldSize: " + QString::number(nField)).toStdString();
            if(nField==nElem){
                typeData = "CELL";
            } else if(nField==nPts){
                typeData = "POINT";
            } else {
                MITK_INFO << "Inconsistent file size";
                break;
            }
            CemrgCommonUtils::AppendScalarFieldToVtk(vtkPath, fi2.baseName(), typeData, field);

            appendScalarFieldReply = QMessageBox::question(NULL, "Question",
                    "Append another scalar field from a file?", QMessageBox::Yes, QMessageBox::No);
        }
    }
}
