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
 * Commandline Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

//Qmitk
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>

//Qt
#include <QFileDialog>
#include <QFileInfo>
#include <QFile>
#include <QCoreApplication>
#include <QDebug>
#include <QDir>
#include <QMessageBox>

//Generic
#include <thread>
#include <chrono>
#include <sys/stat.h>

#include "CemrgCommonUtils.h"
#include "CemrgCommandLine.h"

CemrgCommandLine::CemrgCommandLine() {

    _useDockerContainers = true;
    _dockerimage = "biomedia/mirtk:v1.1.0";

    //Setup panel
    panel = new QTextEdit(0,0);
    QPalette palette = panel->palette();
    palette.setColor(QPalette::Base, Qt::black);
    palette.setColor(QPalette::Text, Qt::red);
    panel->setPalette(palette);
    panel->setReadOnly(true);

    //Setup dialog
    layout = new QVBoxLayout();
    dial = new QDialog(0,0);
    dial->setFixedSize(640, 480);
    dial->setLayout(layout);
    dial->layout()->addWidget(panel);
    dial->show();

    //Setup the process
    process = std::unique_ptr<QProcess>(new QProcess(this));
    process->setProcessChannelMode(QProcess::MergedChannels);
    connect(process.get(), SIGNAL(readyReadStandardOutput()), this, SLOT(UpdateStdText()));
    connect(process.get(), SIGNAL(readyReadStandardError()), this, SLOT(UpdateErrText()));
    connect(process.get(), SIGNAL(finished(int)), this, SLOT(FinishedAlert()));
}

CemrgCommandLine::~CemrgCommandLine() {

    process->close();
    dial->deleteLater();
    panel->deleteLater();
    layout->deleteLater();
}

QDialog* CemrgCommandLine::GetDialog() {

    return dial;
}

/***************************************************************************
 ****************** Execute Plugin Specific Functions **********************
 ***************************************************************************/

QString CemrgCommandLine::ExecuteSurf(QString dir, QString segPath, QString morphOperation, int iter, float th, int blur, int smth) {

    MITK_INFO << "[ATTENTION] SURFACE CREATION: Close -> Surface -> Smooth";

    QString closeOutputPath, surfOutputPath;
    QString outAbsolutePath = "ERROR_IN_PROCESSING";
    closeOutputPath = ExecuteMorphologicalOperation(morphOperation, dir, segPath, "segmentation.s.nii", iter);

    mitk::ProgressBar::GetInstance()->Progress();
    if (QString::compare(closeOutputPath, "ERROR_IN_PROCESSING")!=0) {

        surfOutputPath = ExecuteExtractSurface(dir, closeOutputPath, "segmentation.vtk", th, blur);
        mitk::ProgressBar::GetInstance()->Progress();

        if (QString::compare(surfOutputPath, "ERROR_IN_PROCESSING")!=0) {

            outAbsolutePath = ExecuteSmoothSurface(dir, surfOutputPath, surfOutputPath, smth);
            remove((dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.s.nii").toStdString().c_str());
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            mitk::ProgressBar::GetInstance()->Progress();
        }//_if
    } else {
        mitk::ProgressBar::GetInstance()->Progress(2);
    }//_if

    return outAbsolutePath;
}

QString CemrgCommandLine::ExecuteCreateCGALMesh(QString dir, QString outputName, QString paramsFullPath, QString segmentationName) {

    MITK_INFO << "[ATTENTION] Attempting meshtools3d libraries.";

    QString revertDockerImage = "";
    if (!_dockerimage.contains("meshtools3d", Qt::CaseInsensitive)) {
        MITK_INFO << "Changing docker image name to meshtools3d.";
        revertDockerImage = GetDockerImage(); //get current docker image
        SetDockerImage(QString("alonsojasl/meshtools3d:v1.0"));
    }

    QString executablePath = "";
    QString executableName;
    QDir meshtools3dhome(dir);
    QString outAbsolutePath, outputDirectory, segmentationDirectory;
    QStringList arguments;
    segmentationDirectory = dir + mitk::IOUtil::GetDirectorySeparator();
    outputDirectory = segmentationDirectory + "CGALMeshDir";
    outAbsolutePath = outputDirectory + mitk::IOUtil::GetDirectorySeparator() + outputName;
    outAbsolutePath += ".vtk"; // many outputs are created with meshtools3d. .vtk is the one used in CemrgApp

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(meshtools3dhome.absolutePath());
        arguments << "-f" << meshtools3dhome.relativeFilePath(paramsFullPath);
        arguments << "-seg_dir" << meshtools3dhome.relativeFilePath(segmentationDirectory);;
        arguments << "-seg_name" << segmentationName;
        arguments << "-out_dir" << meshtools3dhome.relativeFilePath(outputDirectory);
        arguments << "-out_name" << outputName;

    } else {

        MITK_INFO << "Using static MESHTOOLS3D libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "M3DLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("M3DLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + "meshtools3d";
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << "-f" << paramsFullPath;
            arguments << "-seg_dir" << segmentationDirectory;;
            arguments << "-seg_name" << segmentationName;
            arguments << "-out_dir" << outputDirectory;
            arguments << "-out_name" << outputName;

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MESHTOOLS3D libraries not found");
            MITK_WARN << "MESHTOOLS3D libraries not found. Please make sure the M3DLib folder is inside the directory:\n\t" + mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    //Setup EnVariable - in windows TBB_NUM_THREADS should be set in the system environment variables
#ifndef _WIN32
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    env.insert("TBB_NUM_THREADS","12");
    process->setProcessEnvironment(env);
#endif

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    //Revert to original docker image (in case the object is used later)
    if (!revertDockerImage.isEmpty())
        SetDockerImage(revertDockerImage);

    if (!successful) {
        if (!_useDockerContainers) {
            MITK_WARN << "MESHTOOLS3D did not produce a good outcome. Trying with the MESHTOOLS3D Docker container.";
            SetUseDockerContainersOn();
            return ExecuteCreateCGALMesh(dir, outputName, paramsFullPath, segmentationName);
        } else {
            MITK_WARN << "MESHTOOLS3D Docker container did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else {
        return outAbsolutePath;
    }
}

void CemrgCommandLine::ExecuteTracking(QString dir, QString imgTimes, QString param, QString output) {

    MITK_INFO << "[ATTENTION] Attempting Registration.";

    QString executablePath = "";
    QString executableName;
    QString commandName = "register";
    QDir mirtkhome(dir);
    QString imgTimesFilePath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    imgTimesFilePath = imgTimes.contains(dir, Qt::CaseSensitive) ? imgTimes : prodPath + imgTimes;
    outAbsolutePath = output.contains(dir, Qt::CaseSensitive) ? output : prodPath + output;

    if (!outAbsolutePath.contains(".dof", Qt::CaseSensitive)) outAbsolutePath += ".dof";

    MITK_INFO << ("[...] IMAGE FILES PATH: " + imgTimesFilePath).toStdString();
    MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << "-images" << mirtkhome.relativeFilePath(imgTimesFilePath);
        if (!param.isEmpty()) arguments << "-parin" << mirtkhome.relativeFilePath(param);
        arguments << "-dofout" << mirtkhome.relativeFilePath(outAbsolutePath);
        arguments << "-threads" << "12";
        arguments << "-verbose" << "3";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << "-images" << imgTimesFilePath;
            if (!param.isEmpty()) arguments << "-parin" << param;
            arguments << "-dofout" << outAbsolutePath;
            arguments << "-threads" << "12";
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            ExecuteTracking(dir, imgTimes, param, output);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

void CemrgCommandLine::ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth) {

    //Time Smoothness
    int fctTime = 10;
    noFrames *= smooth;
    if (smooth == 2) {
        fctTime = 5;
    } else if (smooth == 5) {
        fctTime = 2;
    }

    QString output = dir + mitk::IOUtil::GetDirectorySeparator() + "transformed-";
    QString thisOutput;
    for (int i=0; i<noFrames; i++) {
        thisOutput = output + QString::number(i)+".vtk";
        ExecuteTransformationOnPoints(dir, inputMesh, thisOutput, dofin, iniTime);
        iniTime += fctTime;
        mitk::ProgressBar::GetInstance()->Progress();
    }
}

void CemrgCommandLine::ExecuteRegistration(QString dir, QString fixed, QString moving, QString transformFileName, QString modelname) {

    MITK_INFO << "[ATTENTION] Attempting Registration.";

    //lge : fixed   ||   mra : moving
    QString executablePath = "";
    QString executableName;
    QString commandName = "register";
    QDir mirtkhome(dir);
    QString fixedfullpath, movingfullpath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    fixedfullpath = fixed.contains(dir, Qt::CaseSensitive) ? fixed : prodPath + fixed;
    movingfullpath = moving.contains(dir, Qt::CaseSensitive) ? moving : prodPath + moving;
    outAbsolutePath = transformFileName.contains(dir, Qt::CaseSensitive) ? transformFileName : prodPath + transformFileName;

    if (!fixedfullpath.contains(".nii", Qt::CaseSensitive)) fixedfullpath += ".nii";
    if (!movingfullpath.contains(".nii", Qt::CaseSensitive)) movingfullpath += ".nii";
    if (!outAbsolutePath.contains(".dof", Qt::CaseSensitive)) movingfullpath += ".dof";

    MITK_INFO << ("[...] MOVING (source): " + movingfullpath).toStdString();
    MITK_INFO << ("[...] FIXED (target): " + fixedfullpath).toStdString();
    MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << mirtkhome.relativeFilePath(movingfullpath); // input1
        arguments << mirtkhome.relativeFilePath(fixedfullpath); // input2
        arguments << "-dofout" << mirtkhome.relativeFilePath(outAbsolutePath);
        arguments << "-model" << modelname;
        arguments << "-verbose" << "3";
        arguments << "-color";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << movingfullpath;
            arguments << fixedfullpath;
            arguments << "-dofout" << outAbsolutePath;
            arguments << "-model" << modelname;
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    MITK_INFO << ("Performing a " + modelname + " registration").toStdString();
    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            ExecuteRegistration(dir, fixed, moving, transformFileName, modelname);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

void CemrgCommandLine::ExecuteTransformation(QString dir, QString imgname, QString regname, QString transformFileFullPath) {

    MITK_INFO << "[ATTENTION] Attempting Image Transformation.";

    QString executablePath = "";
    QString executableName;
    QString commandName = "transform-image";
    QDir mirtkhome(dir);
    QString dofpath, imgNamefullpath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    imgNamefullpath = imgname.contains(dir, Qt::CaseSensitive) ? imgname : prodPath + imgname;
    outAbsolutePath = regname.contains(dir, Qt::CaseSensitive) ? regname : prodPath + regname;
    dofpath = transformFileFullPath.contains(dir, Qt::CaseSensitive) ? transformFileFullPath : prodPath + transformFileFullPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + imgNamefullpath).toStdString();
    MITK_INFO << ("[...] INPUT DOF: " + dofpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << mirtkhome.relativeFilePath(imgNamefullpath); //input
        arguments << mirtkhome.relativeFilePath(outAbsolutePath); //output
        arguments << "-dof" << mirtkhome.relativeFilePath(dofpath);
        arguments << "-verbose" << "3";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << imgNamefullpath; //input
            arguments << outAbsolutePath; //output
            arguments << "-dof" << dofpath;
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            ExecuteTransformation(dir, imgname, regname, transformFileFullPath);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

void CemrgCommandLine::ExecuteSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString transformFileName, bool transformThePoints) {

    MITK_INFO << "[ATTENTION] Attempting INIT-DOF.";

    QString aPath, executableName, commandName, sourceMeshPath, targetMeshPath, outAbsolutePath, prodPath;
    QStringList arguments;
    QDir mirtkhome(dir);

    aPath = "";
    commandName = "init-dof"; //simple translation
    sourceMeshPath = sourceMeshP.contains(dir, Qt::CaseSensitive) ? sourceMeshP : prodPath + sourceMeshP;
    targetMeshPath = targetMeshP.contains(dir, Qt::CaseSensitive) ? targetMeshP : prodPath + targetMeshP;
    outAbsolutePath = transformFileName.contains(dir, Qt::CaseSensitive) ? transformFileName : prodPath + transformFileName;

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        aPath = "/usr/local/bin/";
#endif
        executableName = aPath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(),  commandName);
        arguments << mirtkhome.relativeFilePath(outAbsolutePath);
        arguments << "-translations" << "-norotations" << "-noscaling" << "-noshearing";
        if (transformThePoints) {
            arguments << "-displacements";
            arguments << mirtkhome.relativeFilePath(sourceMeshPath);
            arguments << mirtkhome.relativeFilePath(targetMeshPath);
        } else {
            arguments << "-source" << mirtkhome.relativeFilePath(sourceMeshPath);
            arguments << "-target" << mirtkhome.relativeFilePath(targetMeshPath);
        }//_if
        arguments << "-verbose" << "3";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = aPath + mitk::IOUtil::GetDirectorySeparator() + commandName;
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        QDir apathd(aPath);

        if (apathd.exists()) {

            process->setWorkingDirectory(aPath);
            arguments << outAbsolutePath;
            arguments << "-translations" << "-norotations" << "-noscaling" << "-noshearing";
            if (transformThePoints) {
                arguments << "-displacements";
                arguments << sourceMeshPath;
                arguments << targetMeshPath;
            } else {
                arguments << "-source" << sourceMeshPath;
                arguments << "-target" << targetMeshPath;
            }
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            ExecuteSimpleTranslation(dir, sourceMeshP, targetMeshP, transformFileName, transformThePoints);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

/***************************************************************************
 ****************** Execute MIRTK Specific Functions **********************
 ***************************************************************************/

QString CemrgCommandLine::ExecuteMorphologicalOperation(QString operation, QString dir, QString segPath, QString outputPath, int iter) {

    MITK_INFO << "[ATTENTION] Attempting Pointset transformation.";

    QString executablePath = "";
    QString executableName;
    QString commandName;
    QStringList arguments;

    if (QString::compare(operation, "dilate", Qt::CaseInsensitive)==0) {
        commandName = "dilate-image";
    } else if (QString::compare(operation, "erode", Qt::CaseInsensitive)==0) {
        commandName = "erode-image";
    } else if (QString::compare(operation, "open", Qt::CaseInsensitive)==0) {
        commandName = "open-image";
    } else if (QString::compare(operation, "close", Qt::CaseInsensitive)==0) {
        commandName = "close-image";
    } else {
        MITK_ERROR << ("Morphological operation: " + operation + " misspelled or not supported.").toStdString();
        return "ERROR_IN_PROCESSING";
    }

    QDir mirtkhome(dir);
    QString inputImgFullPath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();

    inputImgFullPath = segPath.contains(dir, Qt::CaseSensitive) ? segPath : prodPath + segPath;
    outAbsolutePath = outputPath.contains(dir, Qt::CaseSensitive) ? outputPath : prodPath + outputPath;

    MITK_INFO << ("[...] OPERATION: " + operation).toStdString();
    MITK_INFO << ("[...] INPUT IMAGE: " + inputImgFullPath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {
        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << mirtkhome.relativeFilePath(inputImgFullPath);
        arguments << mirtkhome.relativeFilePath(outAbsolutePath);
        arguments << "-iterations" << QString::number(iter);

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << inputImgFullPath;
            arguments << outAbsolutePath;
            arguments << "-iterations" << QString::number(iter);

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            return ExecuteMorphologicalOperation(operation, dir, segPath, outputPath, iter);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else {
        return outAbsolutePath;
    }
}

QString CemrgCommandLine::ExecuteExtractSurface(QString dir, QString segPath, QString outputPath,float th, int blur) {

    MITK_INFO << "[ATTENTION] Attempting Surface extraction.";

    QString executablePath = "";
    QString executableName;
    QString commandName;
    commandName = "extract-surface";
    QDir mirtkhome(dir);
    QString inputImgFullPath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    inputImgFullPath = segPath.contains(dir, Qt::CaseSensitive) ? segPath : prodPath + segPath;
    outAbsolutePath = outputPath.contains(dir, Qt::CaseSensitive) ? outputPath : prodPath + outputPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + inputImgFullPath).toStdString();
    MITK_INFO << ("[...] OUTPUT MESH: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << mirtkhome.relativeFilePath(inputImgFullPath);
        arguments << mirtkhome.relativeFilePath(outAbsolutePath);
        arguments << "-isovalue" << QString::number(th);
        arguments << "-blur" << QString::number(blur);
        arguments << "-ascii";
        arguments << "-verbose" << "3";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << inputImgFullPath;
            arguments << outAbsolutePath;
            arguments << "-isovalue" << QString::number(th);
            arguments << "-blur" << QString::number(blur);
            arguments << "-ascii";
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            return ExecuteExtractSurface(dir, segPath, outputPath,th, blur);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else {
        return outAbsolutePath;
    }
}

QString CemrgCommandLine::ExecuteSmoothSurface(QString dir, QString segPath, QString outputPath, int smth) {

    MITK_INFO << "[ATTENTION] Attempting Surface extraction.";

    QString executablePath = "";
    QString executableName;
    QString commandName;
    commandName = "smooth-surface";
    QDir mirtkhome(dir);
    QString inputMeshFullPath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    inputMeshFullPath = segPath.contains(dir, Qt::CaseSensitive) ? segPath : prodPath + segPath;
    outAbsolutePath = outputPath.contains(dir, Qt::CaseSensitive) ? outputPath : prodPath + outputPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + inputMeshFullPath).toStdString();
    MITK_INFO << ("[...] OUTPUT MESH: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << mirtkhome.relativeFilePath(inputMeshFullPath);
        arguments << mirtkhome.relativeFilePath(outAbsolutePath);
        arguments << "-iterations" << QString::number(smth);
        arguments << "-verbose" << "3";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << inputMeshFullPath;
            arguments << outAbsolutePath;
            arguments << "-iterations" << QString::number(smth);
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            return ExecuteSmoothSurface(dir, segPath, outputPath, smth);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else {
        return outAbsolutePath;
    }
}

void CemrgCommandLine::ExecuteTransformationOnPoints(QString dir, QString meshFullPath, QString outputMeshFullPath, QString transformFileFullPath, double applyingIniTime) {

    MITK_INFO << "[ATTENTION] Attempting Pointset transformation.";

    QString executablePath = "";
    QString executableName;
    QString commandName = "transform-points";
    QDir mirtkhome(dir);
    QString dofpath, inputMeshFullPath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    inputMeshFullPath = meshFullPath.contains(dir, Qt::CaseSensitive) ? meshFullPath : prodPath + meshFullPath;
    outAbsolutePath = outputMeshFullPath.contains(dir, Qt::CaseSensitive) ? outputMeshFullPath : prodPath + outputMeshFullPath;
    dofpath = transformFileFullPath.contains(dir, Qt::CaseSensitive) ? transformFileFullPath : prodPath + transformFileFullPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + inputMeshFullPath).toStdString();
    MITK_INFO << ("[...] INPUT DOF: " + dofpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << mirtkhome.relativeFilePath(inputMeshFullPath); // input
        arguments << mirtkhome.relativeFilePath(outAbsolutePath); // output
        arguments << "-dofin" << mirtkhome.relativeFilePath(dofpath);
        arguments << "-ascii";
        if (applyingIniTime != -100) {
            // -100 is the default value indicating ExecuteApplying is not being called.
            arguments << "-St";
            arguments << QString::number(applyingIniTime);
        }
        arguments << "-verbose" << "3";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << inputMeshFullPath; // input
            arguments << outAbsolutePath; // output
            arguments << "-dofin" << dofpath;
            arguments << "-ascii";
            if (applyingIniTime != -100) {
                // -100 is the default value indicating ExecuteApplying is not being called.
                arguments << "-St";
                arguments << QString::number(applyingIniTime);
            }
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            ExecuteTransformationOnPoints(dir, meshFullPath, outputMeshFullPath, transformFileFullPath);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

void CemrgCommandLine::ExecuteResamplingOnNifti(QString niiFullPath, QString outputNiiFullPath, int isovalue) {

    MITK_INFO << "[ATTENTION] Attempting Image Transformation.";

    QString executablePath = "";
    QString executableName;
    QString commandName = "resample-image";
    QFileInfo inputnii(niiFullPath);
    QString dir = inputnii.absolutePath();
    QDir mirtkhome(dir);
    QString dofpath, imgNamefullpath, outAbsolutePath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    imgNamefullpath = niiFullPath;
    outAbsolutePath = outputNiiFullPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + imgNamefullpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";

        arguments = GetDockerArguments(mirtkhome.absolutePath(), commandName);
        arguments << mirtkhome.relativeFilePath(imgNamefullpath); //input
        arguments << mirtkhome.relativeFilePath(outAbsolutePath); //output
        arguments << "-isotropic" << QString::number(isovalue);
        arguments << "-interp" << "CSpline";
        arguments << "-verbose" << "3";

    } else {

        MITK_INFO << "Using static MIRTK libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") + mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") + mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);

        if (apathd.exists()) {

            process->setWorkingDirectory(executablePath);
            arguments << imgNamefullpath; //input
            arguments << outAbsolutePath; //output
            arguments << "-isotropic" << QString::number(isovalue);
            arguments << "-interp" << "CSpline";
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);

    if (!successful) {
        if (_useDockerContainers) {
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            SetUseDockerContainersOff();
            ExecuteResamplingOnNifti(niiFullPath, outputNiiFullPath,isovalue);
        } else {
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

/***************************************************************************
 ****************** Execute Docker Specific Functions **********************
 ***************************************************************************/

QString CemrgCommandLine::DockerCemrgNetPrediction(QString mra) {

    MITK_INFO << "[CEMRGNET] Attempting prediction using Docker";

    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif
    QFileInfo finfo(mra);
    QDir cemrgnethome(finfo.absolutePath());
    QString inputfilepath = cemrgnethome.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "test.nii";
    QString tempfilepath = cemrgnethome.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "output.nii";
    QString outputfilepath = cemrgnethome.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "LA-cemrgnet.nii";
    bool test;
    QString res;

    if (QFile::exists(inputfilepath)) {
        MITK_INFO << "[CEMRGNET] File test.nii exists.";
        test = true;
    } else {
        MITK_INFO << "[CEMRGNET] Copying file to test.nii";
        test = QFile::copy(finfo.absoluteFilePath(), inputfilepath);
    }//_if

    if (test) {

        QString inputRelativePath = cemrgnethome.relativeFilePath(inputfilepath);
        process->setWorkingDirectory(cemrgnethome.absolutePath());

        //Setup docker
        QString docker = aPath+"docker";
        QString dockerimage = "orodrazeghi/cemrgnet";
        QStringList arguments;
        arguments << "run" << "--rm";
        arguments << "--volume="+cemrgnethome.absolutePath()+":/data";
        arguments << dockerimage;

        bool debugvar=true;
        if (debugvar) {
            MITK_INFO << "[DEBUG] Input path:";
            MITK_INFO << inputfilepath.toStdString();
            MITK_INFO << "[DEBUG] Docker command to run:";
            MITK_INFO << PrintFullCommand(docker, arguments);
        }

        completion = false;
        process->start(docker, arguments);
        CheckForStartedProcess();
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }

        bool test2 = QFile::rename(tempfilepath, outputfilepath);
        if (test2) {
            MITK_INFO << "[CEMRGNET] Prediction and output creation - successful.";
            res = outputfilepath;
        } else if (IsOutputSuccessful(tempfilepath)) {
            MITK_INFO << "[CEMRGNET] Prediction - successful.";
            res = tempfilepath;
        } else {
            MITK_WARN << "[CEMRGNET] Problem with prediction.";
            res = "";
        }//_if

    } else {
        MITK_WARN << "Copying input file to 'test.nii' was unsuccessful.";
        res = "";
    }//_if

    return res;
}

QString CemrgCommandLine::DockerDicom2Nifti(QString path2dicomfolder) {

    MITK_INFO << "[ATTENTION] Attempting alternative DICOM to NIFTI conversion.";

    QDir dicomhome(path2dicomfolder);
    QString outAbsolutePath = "ERROR_IN_PROCESSING";
    QString executablePath = "";
    QString outPath;
    bool successful = false;

    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        outPath = dicomhome.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "NIIs";
        QString executableName = executablePath+"docker";

        QStringList arguments;
        arguments << "run" << "--rm"  << "--volume="+dicomhome.absolutePath()+":/Data";
        arguments << "orodrazeghi/dicom-converter" << ".";
        arguments << "--gantry" << "--inconsistent";

        successful = ExecuteCommand(executableName, arguments, outPath);

    } else {
        MITK_WARN << "Docker must be running for this feature to be used.";
    }//_if

    if (successful) {
        MITK_INFO << "Conversion successful.";
        outAbsolutePath = outPath;
    } else {
        MITK_WARN << "Error with DICOM2NIFTI Docker container.";
    }

    return outAbsolutePath;
}

QString CemrgCommandLine::DockerSurfaceFromMesh(QString dir, QString meshname, QString outname, QString op, QString outputSuffix){
    // Method equivalent to:  meshtool extract surface
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath = "";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);

    outname += (outputSuffix.at(0)=="_") ? outputSuffix : ("_"+outputSuffix);

    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "extract" << "surface";
    arguments << ("-msh="+meshname);
    arguments << ("-op="+op);
    arguments << ("-surf="+outname);

    QString outPath = home.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + outname + ".surf.vtx";

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (successful) {
        MITK_INFO << "Surface extraction successful.";
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << "Error with MESHTOOL Docker container.";
    }
    return outAbsolutePath;
}

QString CemrgCommandLine::DockerExtractGradient(QString dir, QString meshname, QString idatName, QString odatName, bool elemGrad){
    // Method equivalent to:  meshtool extract surface
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath = "";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);

    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "extract" << "gradient";
    arguments << ("-msh="+meshname);
    arguments << ("-idat="+home.relativeFilePath(idatName));
    arguments << ("-odat="+home.relativeFilePath(odatName));
    if(elemGrad){
        arguments << "-mode=1"; // compute element gradient
    }

    QString outPath = home.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + odatName + ".grad.vec";

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (successful) {
        MITK_INFO << "Gradient extraction successful.";
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << "Error with MESHTOOL Docker container.";
    }
    return outAbsolutePath;
}

QString CemrgCommandLine::OpenCarpDockerLaplaceSolves(QString dir, QString meshName, QString outName, QStringList zeroName, QStringList oneName, QStringList regionLabels){
    SetDockerImage("docker.opencarp.org/opencarp/opencarp:latest");
    QString executablePath = "";
    #if defined(__APPLE__)
            executablePath = "/usr/local/bin/";
    #endif
        QString executableName = executablePath+"docker";
        QString outAbsolutePath = "ERROR_IN_PROCESSING";

        QDir home(dir);
        QString outPath = home.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "fibres_pt_1"+mitk::IOUtil::GetDirectorySeparator()+outName;
        QDir outDir(outPath);
        QString outIgbFile = outPath + mitk::IOUtil::GetDirectorySeparator() + "phie.igb";

        MITK_INFO(outDir.mkpath(outPath)) << "Output directory created.";
        if(!outDir.exists()){
            MITK_INFO << ("Error creating directory: " + outPath).toStdString();
        } else{
            QStringList arguments;
            arguments << "run" << "--rm" << ("--volume="+home.absolutePath()+":/shared:z") << "--workdir=/shared";
            arguments << "docker.opencarp.org/opencarp/opencarp:latest";
            arguments << "openCARP";
            arguments << "-ellip_use_pt" << "0" << "-parab_use_pt" << "0";
            arguments << "-parab_options_file";
            arguments << "/usr/local/lib/python3.6/dist-packages/carputils-0.0.0-py3.6-linux-x86_64.egg/carputils/resources/petsc_options/ilu_cg_opts";
            arguments << "-ellip_options_file";
            arguments << "/usr/local/lib/python3.6/dist-packages/carputils-0.0.0-py3.6-linux-x86_64.egg/carputils/resources/petsc_options/amg_cg_opts";
            arguments << "-simID" << home.relativeFilePath(outPath);
            arguments << "-meshname" << meshName;
            arguments << "-experiment" << "2";
            arguments << "-bidomain" << "1";

            MITK_INFO << "[openCarp] setup arguments computed";

            int numStim = zeroName.size() + oneName.size();
            arguments << "-num_stim" << QString::number(numStim);
            for (int ix = 0; ix < zeroName.size(); ix++) {
                arguments << ("-stimulus[" + QString::number(ix) + "].vtx_file") << home.relativeFilePath(zeroName.at(ix));
                arguments << ("-stimulus[" + QString::number(ix) + "].stimtype") << "3"; // =0
            }
            for (int jx = zeroName.size(); jx < numStim; jx++) {
                arguments << ("-stimulus[" + QString::number(jx) + "].vtx_file") << home.relativeFilePath(oneName.at(jx-zeroName.size()));
                arguments << ("-stimulus[" + QString::number(jx) + "].stimtype") << "2"; // =1
            }

            arguments << ("-stimulus[" + QString::number(numStim-1) + "].start") << "0";
            arguments << ("-stimulus[" + QString::number(numStim-1) + "].duration") << "1";
            arguments << ("-stimulus[" + QString::number(numStim-1) + "].strength") << "1";

            MITK_INFO << "[openCarp] stimulus arguments computed";

            arguments << "-num_gregions" << "1";
            arguments << "-gregion[0].name" << "myo";
            arguments << "-gregion[0].g_il" << "1";
            arguments << "-gregion[0].g_it" << "1";
            arguments << "-gregion[0].g_in" << "1";
            arguments << "-gregion[0].g_el" << "1";
            arguments << "-gregion[0].g_et" << "1";
            arguments << "-gregion[0].g_en" << "1";
            arguments << "-gregion[0].num_IDs"<< QString::number(regionLabels.size());
            for (int ix = 0; ix < regionLabels.size(); ix++) {
                arguments << "-gregion[0].ID[" +QString::number(ix)+ "]"<< regionLabels.at(ix);
            }

            bool successful = ExecuteCommand(executableName, arguments, outIgbFile);

            if (successful) {
                MITK_INFO << "Laplace solves generation successful. Creating .dat file";
                QString outPathFile = dir+mitk::IOUtil::GetDirectorySeparator()+meshName+"_"+outName+"_potential.dat";

                arguments.clear();
                arguments << "run" << "--rm" << ("--volume="+home.absolutePath()+":/shared:z") << "--workdir=/shared";
                arguments << "docker.opencarp.org/opencarp/opencarp:latest";
                arguments << "igbextract" << home.relativeFilePath(outIgbFile) << "-O";
                arguments << home.relativeFilePath(outPathFile) << "-o" << "ascii_1pL";
                successful = ExecuteCommand(executableName, arguments, outPathFile);
                if(successful){
                    outAbsolutePath =  outPathFile;
                }
            } else{
                MITK_WARN << "Error with openCARP LAPLACE SOLVES Docker container.";
            }
        }

        return outAbsolutePath;
}

QString CemrgCommandLine::DockerComputeFibres(QString dir, QString meshname, QString lapapex, QString lapepi, QString laplv, QString laprv, QString type, double a_endo, double a_epi, double b_endo, double b_epi){
    SetDockerImage("alonsojasl/cemrg-fibres:v1.0");
    QString executablePath = "";

#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);

    QString outPath = home.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + meshname+"_fibres.lon";
    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "-m" << meshname;
    arguments << "--type" << type;
    arguments << "-a" << home.relativeFilePath(lapapex);
    arguments << "-e" << home.relativeFilePath(lapepi);
    arguments << "-l" << home.relativeFilePath(laplv);
    arguments << "-r" << home.relativeFilePath(laprv);
    arguments << "--alpha_endo" << QString::number(a_endo);
    arguments << "--alpha_epi" << QString::number(a_epi);
    arguments << "--beta_endo" << QString::number(b_endo);
    arguments << "--beta_epi" << QString::number(b_epi);
    arguments << "-o" << meshname+"_fibres.lon";

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (successful) {
        MITK_INFO << "Fibres generation successful.";
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << "Error with FIBRES Docker container.";
    }
    return outAbsolutePath;
}

QString CemrgCommandLine::ExecuteLaplaceSolves(QString dir, QString meshName, QString outputName, QString paramsFullPath, bool optCarp, bool optVtk){
    MITK_INFO << "[ATTENTION] Attempting meshtools3d/lapsolves libraries.";
    SetDockerImage("alonsojasl/meshtools3d-lapsolves:v1.0");

    QString outAbsolutePath = "ERROR_IN_PROCESSING";
    QString executablePath = "";
    QString executableName;
    QDir home(dir);

    QStringList arguments;
    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        executableName = executablePath+"docker";
        arguments = GetDockerArguments(home.absolutePath());
        arguments << "-f" << home.relativeFilePath(paramsFullPath);
        arguments << "-mesh_dir" << home.relativeFilePath(dir);
        arguments << "-mesh_name" << meshName;
        arguments << "-out_dir" << home.relativeFilePath(dir);
        arguments << "-out_name" << outputName;

    } else{
        MITK_INFO << "Using static MESHTOOLS3D libraries.";
        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "M3DLib";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("M3DLib");
        #endif

        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + "lapsolves";
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);
            arguments << "-f" << paramsFullPath;
            arguments << "-mesh_dir" << dir;
            arguments << "-mesh_name" << meshName;
            arguments << "-out_dir" << dir;
            arguments << "-out_name" << outputName;
        } else{
            QMessageBox::warning(NULL, "Please check the LOG", "MESHTOOLS3D libraries not found");
            MITK_WARN << "MESHTOOLS3D libraries not found. Please make sure the M3DLib folder is inside the directory:\n\t" + mitk::IOUtil::GetProgramPath();
        }
    }
    if(optCarp){ // to read meshtool generated surfaces (default=true)
        arguments << "-carp";
    }
    if(!optVtk){ // allows creation of output VTK file with laplace solves
        arguments << "-novtk";
    }
#ifndef _WIN32
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    env.insert("TBB_NUM_THREADS","12");
    process->setProcessEnvironment(env);
#endif

    // if _lap_apex.dat is generated, we assume _lap_epi, _lap_lv , and _lap_rv are too.
    QString outname = outputName + "_lap_apex_potential.dat";
    QString outPath = home.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + outname;

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (!successful) {
        if (!_useDockerContainers) {
            MITK_WARN << "MESHTOOLS3D did not produce a good outcome. Trying with the MESHTOOLS3D Docker container.";
            SetUseDockerContainersOn();
            MITK_INFO << "[ATTENTION] Creating param file to docker-compatible parameter file.";
            MITK_INFO << "[ATTENTION] Docker implemenation does not allow for vtk output creation.";
            QString dockerparam = CemrgCommonUtils::M3dlibLapSolvesParamFile(dir, "docker-param.par", meshName, dir, true);
            return ExecuteLaplaceSolves(dir, meshName, outputName, dockerparam, optCarp, false);
        } else {
            MITK_WARN << "MESHTOOLS3D Docker container did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else {
        return outAbsolutePath;
    }
}

/***************************************************************************
 *********************** Docker Helper Functions ***************************
 ***************************************************************************/

void CemrgCommandLine::SetUseDockerContainers(bool dockerContainersOnOff) {

    QString onoff = dockerContainersOnOff ? "ON" : "OFF";
    MITK_INFO << ("[...] Setting _useDockerContainers variable to: " + onoff).toStdString();
    _useDockerContainers = dockerContainersOnOff;
}

QStringList CemrgCommandLine::GetDockerArguments(QString volume, QString dockerexe) {

    bool mirtkTest = QString::compare(_dockerimage, "biomedia/mirtk:v1.1.0", Qt::CaseSensitive);
    QStringList argumentList;
    argumentList << "run" << "--rm"  << "--volume="+volume+":/data";
    argumentList << _dockerimage;
    if (mirtkTest == 0)
        argumentList << dockerexe;
    return argumentList;
}

/***************************************************************************
 **************************** Helper Functions *****************************
 ***************************************************************************/

bool CemrgCommandLine::CheckForStartedProcess() {

    //CHECK FOR STARTED PROCESS
    //This function prevents freezing of the app when something goes wrong with the Qt process.
    bool startedProcess = false;
    bool debugvar = false;

    if (debugvar) {
        QStringList errinfo = QProcess::systemEnvironment();
        QString errorInfoString = "";
        for (int ix=0; ix < errinfo.size(); ix++)
            errorInfoString += errinfo.at(ix) + " ";
        MITK_INFO << "SYSTEM ENVIRONMENT:";
        MITK_INFO << errorInfoString.toStdString();
    }

    if (process->waitForStarted()) {

        MITK_INFO << "Starting process";
        startedProcess = true;

    } else {

        completion = true;
        MITK_WARN << "[ATTENTION] Process error!";
        MITK_INFO << "STATE:";
        MITK_INFO << process->state();
        MITK_INFO << "ERROR:";
        MITK_INFO << process->error();

    }//_if

    return startedProcess;
}

void CemrgCommandLine::ExecuteTouch(QString filepath) {

#ifdef _WIN32
    MITK_INFO << "[ATTENTION] touch command only necessary on macOS systems. Step ignored.";
#else
    QString commandName;
    QStringList arguments;
    commandName = "touch"; // touch filepath
    arguments << filepath;

    completion = false;
    process->start(executableName, arguments);

    bool processStarted = CheckForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    MITK_INFO(!processStarted) << "[ATTENTION] TOUCH Process never started.";
#endif
}

bool CemrgCommandLine::IsOutputSuccessful(QString outputFullPath) {

    MITK_INFO << "[ATTENTION] Checking for successful output on path:";
    MITK_INFO << outputFullPath.toStdString();

    QFileInfo finfo(outputFullPath);
    bool fileExists = finfo.exists();
    bool fileSizeTest = false;
    bool res = false;

    MITK_INFO << (fileExists ? "File exists." : "Output file not found.");
    if(fileExists){
        fileSizeTest = finfo.size() > 0;
        if (fileSizeTest){
            MITK_INFO << ("File size: " + QString::number(finfo.size())).toStdString();
            res = true;
        } else{
            MITK_INFO << "File empty. Output unsuccessful";
        }
    }

    return res;
}

std::string CemrgCommandLine::PrintFullCommand(QString command, QStringList arguments) {

    bool debugging = true;
    QString argumentList = "";
    for (int ix=0; ix < arguments.size(); ix++)
        argumentList += arguments.at(ix) + " ";

    if (debugging) {
        QString prodPath = QString::fromStdString(mitk::IOUtil::GetProgramPath());
        MITK_INFO << ("Program path: " + prodPath).toStdString();
        ofstream prodFile1;
        prodFile1.open((prodPath + "dockerDebug.txt").toStdString(), ofstream::out | ofstream::app);
        prodFile1 << (command + " " + argumentList).toStdString() << "\n";
        prodFile1.close();
    }//_if

    return (command + " " + argumentList).toStdString();
}

bool CemrgCommandLine::ExecuteCommand(QString executableName, QStringList arguments, QString outputPath) {

    MITK_INFO << PrintFullCommand(executableName, arguments);

    MITK_INFO << ("[ExecuteCommand] Creating empty file at output:" + outputPath).toStdString();
    ExecuteTouch(outputPath);

    completion = false;
    process->start(executableName, arguments);
    bool successful = false;
    bool processStarted = CheckForStartedProcess();

    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    if (processStarted)
        successful = IsOutputSuccessful(outputPath);
#ifdef _WIN32
    MITK_INFO << "[ATTENTION] command not necessary on Windows systems. Step ignored.";
#else
    if (successful) {
        QFile(outputPath).copy(outputPath + "copy");
        QFile(outputPath).remove();
        QFile(outputPath + "copy").rename(outputPath);
    }//_if
#endif

    return successful;
}

/***************************************************************************
 ************************** Protected Slots ********************************
 ***************************************************************************/

void CemrgCommandLine::UpdateStdText() {

    QByteArray data = process->readAllStandardOutput();
    panel->append(QString(data));
}

void CemrgCommandLine::UpdateErrText() {

    QByteArray data = process->readAllStandardError();
    panel->append(QString(data));
}

void CemrgCommandLine::FinishedAlert() {

    completion = true;
    QString data = process->program() + " Completed!";
    panel->append(data);
}
