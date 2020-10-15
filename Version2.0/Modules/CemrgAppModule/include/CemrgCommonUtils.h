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
 * Simple Common Utilities
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgCommonUtils_h
#define CemrgCommonUtils_h

#include <MitkCemrgAppModuleExports.h>
#include <mitkImage.h>
#include <mitkBoundingObject.h>
#include <mitkDataNode.h>
#include <mitkDataStorage.h>
#include <QString>

class MITKCEMRGAPPMODULE_EXPORT CemrgCommonUtils {

public:

    //Cropping Utils
    static mitk::Image::Pointer CropImage();
    static void SetImageToCut(mitk::Image::Pointer imageToCut);
    static void SetCuttingCube(mitk::BoundingObject::Pointer cuttingCube);
    static void SetImageNode(mitk::DataNode::Pointer imageNode);
    static void SetCuttingNode(mitk::DataNode::Pointer cuttingNode);
    static mitk::DataNode::Pointer GetImageNode();
    static mitk::DataNode::Pointer GetCuttingNode();

    //Sampling Utils
    static mitk::Image::Pointer Downsample(mitk::Image::Pointer image, int factor);
    static mitk::Image::Pointer IsoImageResampleReorient(mitk::Image::Pointer image, bool resample=false, bool reorientToRAI=false);
    static mitk::Image::Pointer IsoImageResampleReorient(QString imPath, bool resample=false, bool reorientToRAI=false);

    // convert2nifti
    static bool ConvertToNifti(mitk::BaseData::Pointer oneNode, QString path2file, bool resample=false, bool reorient=false);

    //Mesh Utils
    static mitk::Surface::Pointer LoadVTKMesh(std::string path);
    static QString M3dlibParamFileGenerator(QString dir, QString filename="param-template.par", QString thicknessCalc="0");
    static void ConvertToCarto(std::string vtkPath);
    static void CalculatePolyDataNormals(vtkSmartPointer<vtkPolyData>& pd, bool celldata=true);

    //Carp Utils
    static void OriginalCoordinates(QString imagePath, QString pointPath, QString outputPath, double scaling=1000);
    static void CalculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath);
    static void RegionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath);

    //Generic
    static mitk::DataNode::Pointer AddToStorage(
            mitk::BaseData* data, std::string nodeName, mitk::DataStorage::Pointer ds);

private:
    //Cropping Utils
    static mitk::Image::Pointer imageToCut;
    static mitk::BoundingObject::Pointer cuttingCube;
    static mitk::DataNode::Pointer imageNode;
    static mitk::DataNode::Pointer cuttingNode;
};

#endif // CemrgCommonUtils_h
