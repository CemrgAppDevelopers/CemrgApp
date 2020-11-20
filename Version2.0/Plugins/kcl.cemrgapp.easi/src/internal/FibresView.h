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


#ifndef FibresView_h
#define FibresView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QMessageBox>
#include <vtkIdList.h>
#include <vtkActor.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkStreamTracer.h>
#include <vtkTubeFilter.h>
#include <vtkPointSource.h>
#include <string>
#include <sstream>
#include <CemrgAtriaClipper.h>
#include <CemrgScarAdvanced.h>
#include <mitkUnstructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

#include "QmitkRenderWindow.h"
#include "mitkCommon.h"
#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "mitkSurface.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"

#include "ui_FibresViewControls.h"
#include "ui_FibresViewUITags.h"
#include "ui_FibresViewUIAngles.h"

/**
  \brief FibresView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class FibresView : public QmitkAbstractView {
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:
  static const std::string VIEW_ID;
  static void SetDirectoryFile(const QString directory, const QString basename);

  void SetNames(); // sets all names for calculations
  void ExtractSurfacesProcess();
  void LaplaceSolvesProcess();
  void LaplaceSolvesOpenCarpProcess();
  void GenerateFibresProcess();

  // helper functions
  QString CreateSurfExtractOp(QStringList vect, QStringList joinings);
  std::vector<double> ReadInScalarField(QString fieldPath);
  std::vector<double> ReadInVectorField(QString fieldPath);
  void CheckTag(QString someTag, QString tagname);
  bool IsPreprocessingDone();
  bool IsExtractSurfDone();
  bool IsLaplaceSolvesDone();
  bool IsFibresGenerationDone();

  void FinishedProcess(QString finishedProcessName, QString xtramessage="");

  // visualiser functions
  void ScalarFieldSelector(int vtkIndex);
  void SetStreamTracerForFibres();

  // inline setters
  void SetDefaultTag(QString tagname);
  inline void SetDefaultTagLVendo(){tag_lvendo = "10";};
  inline void SetDefaultTagLVepi(){tag_lvepi = "1";};
  inline void SetDefaultTagLVbase(){tag_lvbase = "2";};
  inline void SetDefaultTagRVendo(){tag_rvendo = "30";};
  inline void SetDefaultTagRVepi(){tag_rvepi = "5";};
  inline void SetDefaultTagRVbase(){tag_rvbase = "4";};
  inline void SetDefaultTagApex(){tag_apex = "3";};

  inline void SetDefaultAngleAlphaEndo(){alpha_endo = 40;};
  inline void SetDefaultAngleAlphaEpi(){alpha_epi = -50;};
  inline void SetDefaultAngleBetaEndo(){beta_endo = -65;};
  inline void SetDefaultAngleBetaEpi(){beta_epi = 25;};

protected slots:
    // functionality slots
    void Preprocessing(); //btn1
    void ExtractSurfaces(); //btn2
    void LaplaceSolves(); //btn3
    void GenerateFibres(); //btn4

    // ui slots
    void CheckForVtk(); //btnX
    void Visualise(); //btnY
    void VisualiseFibres(); //btnZ

protected:
    virtual void OnSelectionChanged(
        berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;
    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;

    Ui::FibresViewControls m_Controls;
    Ui::FibresViewUITags m_UITags;
    Ui::FibresViewUIAngles m_UIAngles;

private:
    void PreSurf();
    void Visualiser(int vtkIndex=0);

    static QString basename; // normally 'mesh'
    static QString directory; // working directory passed from EASI view

    QString dir, name;
    QString elemFull, ptsFull, shiftPts, cogPts, cavElem;
    QString modPathName, modName;
    int nPts, nElem;
    std::vector<int> meshAttributes;

    QString surfApex, surfBase, surfEpi, surfLVendo, surfRVendo;
    QString lapApex, lapEpi, lapLV, lapRV;
    QString imagePath, dilatedCavPath;
    QString fibresFile, vtkPath;
    QString tag_lvendo, tag_lvepi, tag_lvbase, tag_rvendo, tag_rvepi, tag_rvbase, tag_apex;
    double alpha_endo, alpha_epi, beta_endo, beta_epi;
    bool preprocessDone, extractSurfaceDone, laplaceSolvesDone, fibresDone;
    bool fibresVtkCreated;

    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkActor> meshActor;

    // mitk::UnstructuredGrid::Pointer mitkVtkGrid;
    // vtkSmartPointer<vtkUnstructuredGrid> vtkGrid;
    vtkSmartPointer<vtkUnstructuredGridReader> reader;
    vtkSmartPointer<vtkTubeFilter> tubeFilter;

};

#endif // FibresView_h
