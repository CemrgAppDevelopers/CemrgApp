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

#ifndef CemrgCommandLine_h
#define CemrgCommandLine_h

// Qt
#include <memory>
#include <QProcess>
#include <QTextEdit>
#include <QVBoxLayout>
#include <MitkCemrgAppModuleExports.h>

class MITKCEMRGAPPMODULE_EXPORT CemrgCommandLine : public QObject {

    Q_OBJECT
    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)

public:

    CemrgCommandLine();
    ~CemrgCommandLine();
    QDialog* GetDialog();

    //Execute Plugin Specific Functions
    QString ExecuteSurf(QString dir, QString segPath, QString morphOperation="close", int iter=1, float th=0.5, int blur=0, int smth=10);
    QString ExecuteCreateCGALMesh(QString dir, QString outputName, QString paramsFullPath, QString segmentationName="converted.inr");
    QString ExecuteLaplaceSolves(QString dir, QString meshName, QString outputName, QString paramsFullPath, bool optCarp=true, bool optVtk=false);
    void ExecuteTracking(QString dir, QString imgTimes, QString param, QString output="tsffd.dof");
    void ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth);
    void ExecuteRegistration(QString dir, QString fixed, QString moving, QString transformFileName="rigid.dof", QString modelname="Rigid");
    void ExecuteTransformation(QString dir, QString imgNamefullpath, QString regImgNamefullpath, QString transformFileFullPath="rigid.dof");
    void ExecuteSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString transformFileName="init.dof", bool transformThePoints=true);

    //Execute MIRTK Specific Functions
    QString ExecuteMorphologicalOperation(QString operation, QString dir, QString segPath, QString outputPath = "segmentation.s.nii", int iter=1);
    QString ExecuteExtractSurface(QString dir, QString segPath, QString outputPath = "segmentation.vtk", float th=0.5, int blur=0);
    QString ExecuteSmoothSurface(QString dir, QString segPath, QString outputPath, int smth=10);
    void ExecuteTransformationOnPoints(QString dir, QString meshFullPath, QString outputMeshFullPath, QString transformFileFullPath, double applyingIniTime=-100);
    void ExecuteResamplingOnNifti(QString niiFullPath, QString outputNiiFullPath, int isovalue);

    //Execute Docker Specific Functions
    QString DockerCemrgNetPrediction(QString mra);
    QString DockerDicom2Nifti(QString path2dicomfolder);
    QString DockerSurfaceFromMesh(QString dir, QString meshname, QString outname, QString op, QString outputSuffix); // = meshtool extract surface
    QString DockerComputeFibres(QString dir, QString meshname, QString lapapex, QString lapepi, QString laplv, QString laprv, QString type="biv", double a_endo=40, double a_epi=-50, double b_endo=-65, double b_epi=25);


    //Docker Helper Functions
    void SetUseDockerContainers(bool dockerContainersOnOff);
    inline void SetUseDockerContainersOn() {SetUseDockerContainers(true);};
    inline void SetUseDockerContainersOff() {SetUseDockerContainers(false);};
    inline void SetDockerImage(QString dockerimage) {_dockerimage = dockerimage;};
    inline QString GetDockerImage() {return _dockerimage;};
    QStringList GetDockerArguments(QString volume, QString dockerexe = "");

    //Helper Functions
    bool CheckForStartedProcess();
    void ExecuteTouch(QString filepath);
    bool IsOutputSuccessful(QString outputFullPath);
    std::string PrintFullCommand(QString command, QStringList arguments);
    bool ExecuteCommand(QString executableName, QStringList arguments, QString outputPath);

protected slots:

    void UpdateStdText();
    void UpdateErrText();
    void FinishedAlert();

private:

    //Dial and panel
    QDialog* dial;
    QTextEdit* panel;
    QVBoxLayout* layout;

    //QProcess
    bool completion;
    QString _dockerimage;
    bool _useDockerContainers;
    std::unique_ptr<QProcess> process;
};

#endif // CemrgCommandLine_h
