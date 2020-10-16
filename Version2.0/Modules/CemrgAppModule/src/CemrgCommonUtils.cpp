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

//ITK
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkOrientImageFilter.h>
#include <itkImage.h>

//VTK
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkCell.h>

//Qmitk
#include <mitkBoundingObjectCutter.h>
#include <mitkProgressBar.h>
#include <mitkDataNode.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkIOUtil.h>
#include <mitkDataStorage.h>
#include <mitkRenderingManager.h>

//Qt
#include <QMessageBox>
#include <QString>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include "CemrgCommonUtils.h"


mitk::DataNode::Pointer CemrgCommonUtils::imageNode;
mitk::DataNode::Pointer CemrgCommonUtils::cuttingNode;
mitk::Image::Pointer CemrgCommonUtils::imageToCut;
mitk::BoundingObject::Pointer CemrgCommonUtils::cuttingCube;

mitk::Image::Pointer CemrgCommonUtils::CropImage() {

    //Test input objects
    if (imageToCut.IsNull() || cuttingCube.IsNull())
        return NULL;

    //Prepare the cutter
    mitk::BoundingObjectCutter::Pointer cutter = mitk::BoundingObjectCutter::New();
    cutter->SetBoundingObject(cuttingCube);
    cutter->SetInput(imageToCut);
    cutter->AutoOutsideValueOff();

    //Actual cutting
    try {
        cutter->Update();
        mitk::ProgressBar::GetInstance()->Progress();
    } catch (const itk::ExceptionObject& e) {
        std::string message = std::string("The Cropping filter could not process because of: \n ") + e.GetDescription();
        QMessageBox::warning(
                    NULL, "Cropping not possible!", message.c_str(),
                    QMessageBox::Ok, QMessageBox::NoButton, QMessageBox::NoButton);
        return NULL;
    }//try

    //Cutting successful
    mitk::Image::Pointer resultImage = cutter->GetOutput();
    resultImage->DisconnectPipeline();
    resultImage->SetPropertyList(imageToCut->GetPropertyList()->Clone());
    mitk::ProgressBar::GetInstance()->Progress();

    return resultImage;
}

void CemrgCommonUtils::SetImageToCut(mitk::Image::Pointer imageToCut) {

    CemrgCommonUtils::imageToCut = imageToCut;
}

void CemrgCommonUtils::SetCuttingCube(mitk::BoundingObject::Pointer cuttingCube) {

    CemrgCommonUtils::cuttingCube = cuttingCube;
}

void CemrgCommonUtils::SetImageNode(mitk::DataNode::Pointer imageNode) {

    CemrgCommonUtils::imageNode = imageNode;
}

void CemrgCommonUtils::SetCuttingNode(mitk::DataNode::Pointer cuttingNode) {

    CemrgCommonUtils::cuttingNode = cuttingNode;
}

mitk::DataNode::Pointer CemrgCommonUtils::GetImageNode() {

    return imageNode;
}

mitk::DataNode::Pointer CemrgCommonUtils::GetCuttingNode() {

    return cuttingNode;
}

mitk::Image::Pointer CemrgCommonUtils::Downsample(mitk::Image::Pointer image, int factor) {

    typedef itk::Image<short,3> ImageType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NearestInterpolatorType;

    //Cast to ITK
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(image, itkImage);

    ResampleImageFilterType::Pointer downsampler = ResampleImageFilterType::New();
    downsampler->SetInput(itkImage);

    NearestInterpolatorType::Pointer interpolator = NearestInterpolatorType::New();
    downsampler->SetInterpolator(interpolator);

    downsampler->SetDefaultPixelValue(0);

    ResampleImageFilterType::SpacingType spacing = itkImage->GetSpacing();
    spacing *= (double) factor;
    downsampler->SetOutputSpacing(spacing);

    downsampler->SetOutputOrigin(itkImage->GetOrigin());
    downsampler->SetOutputDirection(itkImage->GetDirection());

    ResampleImageFilterType::SizeType size = itkImage->GetLargestPossibleRegion().GetSize();
    for (int i=0; i<3; ++i)
        size[i] /= factor;

    downsampler->SetSize(size);
    downsampler->UpdateLargestPossibleRegion();

    //Save downsampled image
    image = mitk::ImportItkImage(downsampler->GetOutput())->Clone();
    mitk::ProgressBar::GetInstance()->Progress();
    return image;
}

mitk::Image::Pointer CemrgCommonUtils::IsoImageResampleReorient(mitk::Image::Pointer image, bool resample, bool reorientToRAI){
    MITK_INFO(resample) << "Resampling image to be isometric.";
    MITK_INFO(reorientToRAI) << "Doing a reorientation to RAI.";

    typedef itk::Image<short,3> ImageType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
    typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> BSplineInterpolatorType;
    ImageType::Pointer itkInputImage = ImageType::New();
    ImageType::Pointer resampleOutput = ImageType::New();
    ImageType::Pointer outputImage = ImageType::New();
    mitk::CastToItkImage(image, itkInputImage);

    if(resample){
        ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
        BSplineInterpolatorType::Pointer binterp = BSplineInterpolatorType::New();
        binterp->SetSplineOrder(3);
        resampler->SetInterpolator(binterp);

        resampler->SetInput(itkInputImage);
        resampler->SetOutputOrigin(itkInputImage->GetOrigin());
        ImageType::SizeType input_size = itkInputImage->GetLargestPossibleRegion().GetSize();
        ImageType::SpacingType input_spacing = itkInputImage->GetSpacing();
        ImageType::SizeType output_size;
        ImageType::SpacingType output_spacing;
        output_size[0] = input_size[0] * (input_spacing[0] / 1.0);
        output_size[1] = input_size[1] * (input_spacing[1] / 1.0);
        output_size[2] = input_size[2] * (input_spacing[2] / 1.0);
        output_spacing [0] = 1.0;
        output_spacing [1] = 1.0;
        output_spacing [2] = 1.0;
        resampler->SetSize(output_size);
        resampler->SetOutputSpacing(output_spacing);
        resampler->SetOutputDirection(itkInputImage->GetDirection());
        resampler->UpdateLargestPossibleRegion();

        resampleOutput = resampler->GetOutput();
    } else{
        resampleOutput = itkInputImage;
    }

    if(reorientToRAI){
        typedef itk::OrientImageFilter<ImageType,ImageType> OrientImageFilterType;
        OrientImageFilterType::Pointer orienter = OrientImageFilterType::New();
        orienter->UseImageDirectionOn();
        orienter->SetDesiredCoordinateOrientationToAxial(); // RAI
        orienter->SetInput(resampleOutput);
        orienter->Update();
        outputImage = orienter->GetOutput();
    } else{
        outputImage = resampleOutput;
    }

    image = mitk::ImportItkImage(outputImage)->Clone();
    mitk::ProgressBar::GetInstance()->Progress();
    return image;
}

mitk::Image::Pointer CemrgCommonUtils::IsoImageResampleReorient(QString imPath, bool resample,  bool reorientToRAI){
    return CemrgCommonUtils::IsoImageResampleReorient(mitk::IOUtil::Load<mitk::Image>(imPath.toStdString()), resample, reorientToRAI);
};

bool CemrgCommonUtils::ConvertToNifti(mitk::BaseData::Pointer oneNode, QString path2file, bool resample, bool reorient){
    bool successful = false;
    if (oneNode) {
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(oneNode.GetPointer());
        if (image) { //Test if this data item is an image
            image = CemrgCommonUtils::IsoImageResampleReorient(image, resample, reorient);
            mitk::IOUtil::Save(image, path2file.toStdString());
            successful = true;
        } else{
            MITK_INFO << "[...] Problem casting node data to image";
        }
    } else{
        MITK_INFO << "[...] Problem with node";
    }
    return successful;
}


mitk::Surface::Pointer CemrgCommonUtils::LoadVTKMesh(std::string path) {

    try {
        //Load the mesh
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path);
        vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

        //Prepare points for MITK visualisation
        double Xmin = 0, Xmax = 0, Ymin = 0, Ymax = 0, Zmin = 0, Zmax = 0;
        for (int i=0; i<pd->GetNumberOfPoints(); i++) {
            double* point = pd->GetPoint(i);
            point[0] = -point[0];
            point[1] = -point[1];
            pd->GetPoints()->SetPoint(i, point);
            //Find mins and maxs
            if (i==0) {
                Xmin = point[0];
                Xmax = point[0];
                Ymin = point[1];
                Ymax = point[1];
                Zmin = point[2];
                Zmax = point[2];
            } else {
                if (point[0]<Xmin) Xmin = point[0];
                if (point[0]>Xmax) Xmax = point[0];
                if (point[1]<Ymin) Ymin = point[1];
                if (point[1]>Ymax) Ymax = point[1];
                if (point[2]<Zmin) Zmin = point[2];
                if (point[2]>Zmax) Zmax = point[2];
            }//_if
        }//_for
        double bounds[6] = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
        surface->GetGeometry()->SetBounds(bounds);

        return surface;

    } catch (...) {
        return mitk::Surface::New();
    }//_catch
}

QString CemrgCommonUtils::M3dlibParamFileGenerator(QString dir, QString filename, QString thicknessCalc) {

    QString path2file = dir + mitk::IOUtil::GetDirectorySeparator() + filename;
    QFile fi(path2file);

    if (thicknessCalc.compare("0", Qt::CaseSensitive)!=0 && thicknessCalc.compare("1", Qt::CaseSensitive)!=0) {
        MITK_INFO << "Thickness calculation set to default (OFF)";
        thicknessCalc = "0";
    }

    if (fi.open(QFile::WriteOnly | QFile::Truncate)) {
        QTextStream out(&fi);
        out << "[segmentation]" << "\n\n";
        out << "seg_dir" << "=" << "./example" << "\n";
        out << "seg_name" << "=" << "converted.inr" << "\n";
        out << "mesh_from_segmentation" << "=" << "1" << "\n\n";

        out << "[meshing]" << "\n\n";
        out << "readTheMesh" << "=" << "0" << "\n";
        out << "mesh_dir" << "=" << "." << "\n";
        out << "mesh_name" << "=" << "mesh" << "\n\n";

        out << "facet_angle" << "=" << "30"<< "\n";
        out << "facet_size" << "=" << "5.0" << "\n";
        out << "facet_distance" << "=" << "4"<< "\n";
        out << "cell_rad_edge_ratio" << "=" << "2.0" << "\n";
        out << "cell_size" << "=" << "1.0" << "\n\n";

        out << "rescaleFactor" << "=" << "1.0  # rescaling for carp and vtk output" << "\n\n";

        out << "[laplacesolver]" << "\n\n";
        out << "abs_toll" << "=" << "1e-6 # Also for evaluating the thickness" << "\n";
        out << "rel_toll" << "=" << "1e-6" << "\n";
        out << "itr_max" << "=" << "500" << "\n";
        out << "dimKrilovSp" << "=" << "150" << "\n";
        out << "verbose" << "=" << "0" << "\n\n";

        out << "[output]" << "\n\n";
        out << "outdir" << "=" << "." << "\n";
        out << "name" << "=" << "imgmesh" << "\n\n";

        out << "out_medit" << "=" << "0" << "\n";
        out << "out_carp" << "=" << "1" << "\n";
        out << "out_carp_binary" << "=" << "0" << "\n";
        out << "out_vtk" << "=" << "1" << "\n";
        out << "out_vtk_binary" << "=" << "0" << "\n";
        out << "out_potential" << "=" << "0" << "\n";
        out << "debug_output" << "=" << "0" << "\n";
        out << "debug_frequency" << "=" << "10000" << "\n\n";

        out << "[others]" << "\n\n";
        out << "eval_thickness" << "=" << thicknessCalc << "\n";
        out << "thickalgo" << "=" << "1" << "\n"; //#1: Martin Bishop Algorithm; 2: Cesare Corrado Algorithm
        out << "swapregions" << "=" << "1" << "\n";
        out << "verbose" << "=" << "0" << "\n";

        return path2file;

    } else {
        MITK_WARN << ("File " + path2file + "not created.").toStdString();
        return "ERROR_IN_PROCESSING";
    }
}

void CemrgCommonUtils::ConvertToCarto(std::string vtkPath) {

    //Read vtk from the file
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(vtkPath.c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();

    //Output path
    std::string outputPath = vtkPath.substr(0, vtkPath.size()-4);
    outputPath = outputPath + "-carto.vtk";

    //File
    ofstream cartoFile;
    cartoFile.open(outputPath);

    //Header
    cartoFile << "# vtk DataFile Version 3.0\n";
    cartoFile << "PatientData Anon Anon 00000000\n";
    cartoFile << "ASCII\n";
    cartoFile << "DATASET POLYDATA\n";

    //Points
    cartoFile << "POINTS\t" << pd->GetNumberOfPoints() << "\tfloat\n";
    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        cartoFile << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    cartoFile << "\n";

    //Cells
    cartoFile << "POLYGONS\t";
    cartoFile << pd->GetNumberOfCells() << "\t";
    cartoFile << pd->GetNumberOfCells()*4 << "\n";
    for (int i=0; i<pd->GetNumberOfCells(); i++) {
        vtkCell* cell = pd->GetCell(i);
        vtkIdList* list = cell->GetPointIds();
        cartoFile << "3";
        for (int j=0; j<list->GetNumberOfIds(); j++)
            cartoFile << " " << list->GetId(j);
        cartoFile << "\n";
    }

    //Point data
    vtkSmartPointer<vtkFloatArray> pointData = vtkSmartPointer<vtkFloatArray>::New();
    try {
        pointData = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetScalars());
    } catch (...) {
        MITK_WARN << "Storing point data failed! Check your input";
        return;
    }//_try

    float min = pointData->GetRange()[0];
    float max = pointData->GetRange()[1];

    MITK_INFO << "Storing point data, number of tuples: " << pointData->GetNumberOfTuples();
    MITK_INFO << "Storing point data, number of components: " << pointData->GetNumberOfComponents();

    if (pointData->GetNumberOfTuples() != 0) {

        cartoFile << "\nPOINT_DATA\t";
        cartoFile << pointData->GetNumberOfTuples() << "\n";

        if (pointData->GetNumberOfComponents() == 1) {

            cartoFile << "SCALARS scalars float\n";
            cartoFile << "LOOKUP_TABLE default\n";
            for (int i=0; i<pointData->GetNumberOfTuples(); i++) {
                double value = static_cast<double>(pointData->GetTuple1(i));
                value = (value - min) / (max - min);
                std::stringstream stream;
                stream << std::fixed << std::setprecision(2) << value;
                cartoFile << stream.str() << "\n";
            }
            cartoFile << "\n";

        } else {

            for (int i=0; pointData->GetNumberOfComponents(); i++) {

                cartoFile << "SCALARS " << "scalars" << i << " float\n";
                cartoFile << "LOOKUP_TABLE default\n";
                for (int j=0; j<pointData->GetNumberOfTuples(); j++)
                    cartoFile << pointData->GetTuple(j)[i] << " ";
                cartoFile << "\n";

            }
        }
    }//_point_data

    MITK_INFO << "Storing lookup table, min/max scalar values: " << min << " " << max;

    cartoFile << "LOOKUP_TABLE lookup_table 3\n";
    cartoFile << "0.0 0.0 1.0 1.0" << "\n";
    cartoFile << "1.0 1.0 1.0 1.0" << "\n";
    cartoFile << "0.0 1.0 0.0 1.0" << "\n";
    cartoFile.close();
}

void CemrgCommonUtils::CalculatePolyDataNormals(vtkSmartPointer<vtkPolyData>& pd, bool celldata){
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    vtkSmartPointer<vtkPolyData> tempPD = vtkSmartPointer<vtkPolyData>::New();
    tempPD->DeepCopy(pd);

    if(celldata){
        normals->ComputeCellNormalsOn();
    } else{ // pointdata
        normals->ComputePointNormalsOn();
    }

    normals->SetInputData(tempPD);
    normals->SplittingOff();
    normals->Update();
    pd = normals->GetOutput();
}

mitk::DataNode::Pointer CemrgCommonUtils::AddToStorage(
        mitk::BaseData* data, std::string nodeName, mitk::DataStorage::Pointer ds) {

    if (!data) return mitk::DataNode::New();

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(nodeName);
    ds->Add(node);

    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(ds);
    return node;
}

//Carp Utils - operations with .elem and .pts files
void CemrgCommonUtils::OriginalCoordinates(QString imagePath, QString pointPath, QString outputPath, double scaling){
    if(QFileInfo::exists(imagePath) && QFileInfo::exists(pointPath)){
        typedef itk::Image<uint8_t,3> ImageType;
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(imagePath.toStdString());
        ImageType::Pointer itkInput = ImageType::New();
        ImageType::PointType origin;
        mitk::CastToItkImage(image, itkInput);
        origin = itkInput->GetOrigin();

        std::ifstream pointFileRead;

        int nPts;
        pointFileRead.open(pointPath.toStdString());
        pointFileRead >> nPts;

        double x,y,z;

        std::ofstream outputFileWrite;
        outputFileWrite.open(outputPath.toStdString());
        outputFileWrite << nPts << std::endl;

        double xt=0.0, yt=0.0, zt=0.0;
        for(int iPt=0; iPt < nPts; iPt++){
            pointFileRead >> x;
            pointFileRead >> y;
            pointFileRead >> z;
            if (pointFileRead.eof()) {
                MITK_WARN << " WARNING!: File ended prematurely " << endl;
                break;
            }

            xt = x+(origin[0]*scaling);
            yt = y+(origin[1]*scaling);
            zt = z+(origin[2]*scaling);

            outputFileWrite << std::fixed << xt << " ";
            outputFileWrite << std::fixed << yt << " ";
            outputFileWrite << std::fixed << zt << std::endl;
        }
        pointFileRead.close();
        outputFileWrite.close();
        MITK_INFO << ("Saved to file: " + outputPath).toStdString();

    } else{
        MITK_ERROR(QFileInfo::exists(imagePath)) << ("Could not read file" + imagePath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("Could not read file" + pointPath).toStdString();
    }

}

void CemrgCommonUtils::CalculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath){
    if(QFileInfo::exists(elemPath) && QFileInfo::exists(pointPath)){
        FILE* pointFileRead = fopen(pointPath.toStdString().c_str(), "r");
        FILE* elemFileRead = fopen(elemPath.toStdString().c_str(), "r");
        int nPts, nElem;

        fscanf(pointFileRead, "%d\n", &nPts);

        double* pts_array = (double*)malloc(nPts * 3 * sizeof(double));
    	if (pts_array == NULL) {
    		MITK_ERROR << "pts_array malloc FAIL";
    		return;
    	}

        MITK_INFO << "Beginning input .pts file";
    	for (int i = 0; i < nPts; i++) {
    		double* loc = pts_array + 3*i;
    		fscanf(pointFileRead, "%lf %lf %lf\n", loc, loc + 1, loc + 2);
    	}
    	MITK_INFO << "Completed input .pts file";
    	fclose(pointFileRead);

        fscanf(elemFileRead, "%d\n", &nElem);

        std::ofstream outputFileWrite;
        outputFileWrite.open(outputPath.toStdString());
        outputFileWrite << nElem << " 3"<< std::endl;

        MITK_INFO << "Beginning input .elem file, simultaneous output";
    	for (int i = 0; i < nElem; i++) {
    // 		int* loc = elem_array + 4*i;
    		int p1, p2, p3, p4, region;
    		fscanf(elemFileRead, "Tt %d %d %d %d %d\n", &p1, &p2, &p3, &p4, &region);

    		// Calculate and output cog
    		double x = 0.0, y = 0.0, z = 0.0;

    		double *loc = pts_array + 3*p1;
    		x += loc[0];
    		y += loc[1];
    		z += loc[2];

    		loc = pts_array + 3*p2;
    		x += loc[0];
    		y += loc[1];
    		z += loc[2];

    		loc = pts_array + 3*p3;
    		x += loc[0];
    		y += loc[1];
    		z += loc[2];

    		loc = pts_array + 3*p4;
    		x += loc[0];
    		y += loc[1];
    		z += loc[2];

    		x /= 4.0 * 1000;
    		y /= 4.0 * 1000;
    		z /= 4.0 * 1000;

            outputFileWrite << x << std::endl;
            outputFileWrite << y << std::endl;
            outputFileWrite << z << std::endl;

    	}
    	MITK_INFO << "Completed input .elem file";

    	fclose(elemFileRead);
        outputFileWrite.close();
        MITK_INFO << "Completed input .elem file";


    } else{
        MITK_ERROR(QFileInfo::exists(elemPath)) << ("Could not read file" + elemPath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("Could not read file" + pointPath).toStdString();
    }
}

void CemrgCommonUtils::RegionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath){
    if(QFileInfo::exists(bpPath) && QFileInfo::exists(pointPath) && QFileInfo::exists(elemPath)){
        typedef itk::Image<uint8_t,3> ImageType;
        typedef itk::Index<3> IndexType;
        // ScarImage image(bpPath.toStdString());
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(bpPath.toStdString());
        ImageType::Pointer itkInput = ImageType::New();
        ImageType::SpacingType spacing;
        ImageType::RegionType region;
        ImageType::PointType origin;
        ImageType::SizeType size;

        mitk::CastToItkImage(image, itkInput);
        origin = itkInput->GetOrigin();
        spacing = itkInput->GetSpacing();
        region = itkInput->GetLargestPossibleRegion();
        size = region.GetSize();

        double Tx[4][4];
        int max_i = size[0]-1;
        int max_j = size[1]-1;
        int max_k = size[2]-1;
        int min_i = 0;
        int min_j = 0;
        int min_k = 0;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                Tx[i][j] = (i==j) ? 1.0 : 0.0;
            }
        }

        std::ofstream outputFileWrite;
        std::ifstream cogFileRead, elemFileRead;
        cogFileRead.open(pointPath.toStdString());

        int nElemCOG, dim, count;
        double x, y, z;

        cogFileRead >> nElemCOG;
        cogFileRead >> dim;


        outputFileWrite.open(outputPath.toStdString());
        outputFileWrite << nElemCOG << std::endl;

        MITK_INFO << ("Number of elements (COG file):" + QString::number(nElemCOG)).toStdString();
        MITK_INFO << ("Dimension= " + QString::number(dim)).toStdString();

        elemFileRead.open(elemPath.toStdString());
        int nElem;

        elemFileRead >> nElem;
        if(nElem != nElemCOG){
            MITK_ERROR << "Number of elements in files are not consistent.";
        }

        count = 0;
        char type[2];
        int nodes[4], imregion, newRegion;
        int newRegionCount = 0;

        for(int iElem=0; iElem < nElemCOG; iElem++){
            cogFileRead >> x;
            cogFileRead >> y;
            cogFileRead >> z;
            if(cogFileRead.eof()){
                MITK_WARN << "File ended prematurely";
                break;
            }

            elemFileRead >> type;
            elemFileRead >> nodes[0];
            elemFileRead >> nodes[1];
            elemFileRead >> nodes[2];
            elemFileRead >> nodes[3];
            elemFileRead >> imregion;

            // checking point belonging to imregion (cm2carp/carp_scar_map::inScar())
            double xt, yt, zt;
            xt = Tx[0][0] * x + Tx[0][1] * y + Tx[0][2] * z + Tx[0][3]- origin[0] + spacing[0]/2;
            yt = Tx[1][0] * x + Tx[1][1] * y + Tx[1][2] * z + Tx[1][3]- origin[1] + spacing[0]/2;
            zt = Tx[2][0] * x + Tx[2][1] * y + Tx[2][2] * z + Tx[2][3]- origin[2] + spacing[0]/2;

            int i,j,k;
            i = static_cast<int>(xt/spacing[0]);
            j = static_cast<int>(yt/spacing[1]);
            k = static_cast<int>(zt/spacing[2]);

            if ( i<0 || j<0 || k<0 ) {
                min_i = std::min(i, min_i);
                min_j = std::min(j, min_j);
                min_k = std::min(k, min_k);
                newRegion = 0;
            }
            else if ( (unsigned) i > size[0]-1 || (unsigned) j > size[1]-1 || (unsigned) k > size[2]-1 ) {
                max_i = std::max(i, max_i);
                max_j = std::max(j, max_j);
                max_k = std::max(k, max_k);
                newRegion = 0;

            } else {
                IndexType index = {{i, j, k}};
                newRegion = itkInput->GetPixel(index);
            }

            if(newRegion != 0){
                imregion = newRegion;
                newRegionCount++;
            }

            outputFileWrite << type << " ";
            outputFileWrite << nodes[0] << " ";
            outputFileWrite << nodes[1] << " ";
            outputFileWrite << nodes[2] << " ";
            outputFileWrite << nodes[3] << " ";
            outputFileWrite << imregion << std::endl;

            count++;
        }

        outputFileWrite.close();

        MITK_INFO << ("Number of element COG read: " + QString::number(count)).toStdString();
        MITK_INFO << ("Number of new regions determined: " + QString::number(newRegionCount)).toStdString();

        if ( min_i < 0 || min_j < 0 || min_k < 0 ) {
            MITK_WARN << "WARNING: The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region.";
            MITK_WARN << "If scar does lie in this region, then you need to pad the image at the start by (in pixel space):";
            MITK_WARN << (QString::number(-min_i) + " " + QString::number(-min_j) + " " + QString::number(-min_k)).toStdString();
            MITK_WARN << "And add the following transformation to the TransMatFile (in geometric space):";
            MITK_WARN << (QString::number(-min_i*spacing[0]) + " " + QString::number(-min_j*spacing[1]) + " " + QString::number(-min_k*spacing[2])).toStdString();
        }
        if ( max_i > int(size[0]-1) || max_j > int(size[1]-1) || max_k > int(size[2]-1) ) {
            MITK_WARN << "WARNING: The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region.";
            MITK_WARN << "If scar does lie in this region, then you need to pad the image at the end by (in pixel space):";
            MITK_WARN << (QString::number(max_i-(size[0]-1)) + " " + QString::number(max_j-(size[1]-1)) + " " + QString::number(max_k-(size[2]-1))).toStdString();
            MITK_WARN << "No need to change TransMatFile";
        }

    } else{
        MITK_ERROR(QFileInfo::exists(bpPath)) << ("File does not exist: " + bpPath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("File does not exist: " + pointPath).toStdString();
    }
}
