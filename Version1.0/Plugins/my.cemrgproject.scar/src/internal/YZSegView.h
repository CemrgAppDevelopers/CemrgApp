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
 * Scar Segmentation Plugin for MITK (Created in July 2016)
 *
 * It has a few different algotithms to choose from:
 *
 *     Full Width Half Maximum (FWHM) method,
 *     4, 6 and customised standard deviation methods
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * yan.zhao@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef YZSegView_h
#define YZSegView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <mitkSurface.h>
#include <CemrgScar3D.h>
#include "ui_YZSegViewControls.h"

/*!
  \brief YZSegView
  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.
  \sa QmitkFunctionality
  \ingroup ${plugin_target}_internal
*/
class YZSegView : public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;

    virtual void CreateQtPartControl(QWidget *parent);

protected slots:

    /// \brief Called when the user clicks the GUI button
    void LoadDICOM();
    void ProcessIMGS();
    void ConvertNII();
    void ResampIMGS();
    void SegmentIMGS();
    void SaveSEG();
    void ScarSeg();
    void ScarSeg_FWHM();
    void ScarSeg_SD();
    void ScarSeg_4SD();
    void ScarSeg_6SD();
    void ScarSeg_CustomisedSD();
    void ScarSeg_save();

protected:

    virtual void SetFocus();
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes);

    Ui::YZSegViewControls m_Controls;

private:

    QString fileName;
    QString directory;
    std::unique_ptr<CemrgScar3D> scar;
};

#endif // YZSegView_h
