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

#ifndef QmitkCemrgAppCommonTools_h
#define QmitkCemrgAppCommonTools_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include "ui_QmitkCemrgAppCartoExport.h"
#include "ui_QmitkCemrgAppCommonToolsControls.h"

/**
  \brief QmitkCemrgAppCommonTools
  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.
  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class QmitkCemrgAppCommonTools : public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;

protected slots:

    /// \brief Called when the user clicks the GUI button
    void LoadMesh();
    void ConvertToCarto();
    void ConvertToCartoUIUpdate();
    void ConvertToCartoUITextUpdate();
    void ConvertCarpToVtk();

protected:

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    virtual void OnSelectionChanged(
            berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::QmitkCemrgAppCommonToolsControls m_Controls;
    Ui::QmitkCemrgAppCartoExport m_CartoUIThresholding;
};

#endif // QmitkCemrgAppCommonTools_h
