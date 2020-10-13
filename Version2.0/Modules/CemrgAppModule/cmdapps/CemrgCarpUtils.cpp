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
CEMRG CMD CARP UTILS
This app serves as a template for the command line apps to be implemented
in the framework.
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkSurface.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkCommandLineParser.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkImplicitBoolean.h>
#include <vtkImplicitVolume.h>
#include <vtkImplicitDataSet.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkIdList.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkDataSetSurfaceFilter.h>

// ITK
#include <itkPoint.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkOrientImageFilter.h>

// Qt
#include <QtDebug>
#include <QString>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>
#include <numeric>

#include <CemrgScar3D.h>
#include <CemrgCommandLine.h>

#include <algorithm>
#include <string>
// #include <math.h>

typedef itk::Image<uint8_t,3> ImageType;
typedef itk::Index<3> IndexType;
typedef itk::ImageRegionIterator<ImageType> IteratorType;

int checkAlgorithm(std::string algo);
bool checkAlgorithmParameters(int algorithmNumber, std::string inImage, std::string inPts, std::string inElem, std::string inBp);

// main functions
bool originalCoordinates(QString imagePath, QString pointPath, QString outputPath, bool v);
bool calculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath, bool v);
bool regionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath, bool v);

// helper functions
int inRegion(ImageType::Pointer image, double x, double y, double z);
bool checkfile(QString filepath);

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("EASI processing");
    parser.setTitle("CARP utils for CemrgApp");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("CARP utils for CemrgApp Command-line App.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::InputFile,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
                "method", "m", mitkCommandLineParser::String,
                "Method to run.", "Methods include: original_coords (coords), calc_cog (cog) and region_map (region).",
                us::Any(), false);
    parser.addArgument( // optional
                "image", "i", mitkCommandLineParser::InputFile,
                "Input image file", "Full path to image file (.nii). Used on 'coords' and 'region' programs.");
    parser.addArgument( // optional
                "points", "pts", mitkCommandLineParser::InputFile,
                "Input points file", "Full path to points file (.pts). Used on all programs.");
    parser.addArgument( // optional
                "elements", "elem", mitkCommandLineParser::InputFile,
                "Input elements file", "Full path to elements file (.elem). Used on 'cog' and 'region' programs.");
    parser.addArgument(// optional
                "bloodpool", "bp", mitkCommandLineParser::String,
                "Input image of dilated bloodpool", "Full path to image file (.nii). Used on 'region' programs.");
    parser.addArgument(// optional
                "output", "o", mitkCommandLineParser::String,
                "Output file name", "Where to save the output (default: output.nii).");
    parser.addArgument( // optional
                "verbose", "v", mitkCommandLineParser::Bool,
                "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["method"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    auto algo = us::any_cast<std::string>(parsedArgs["method"]);
    int algorithmNumber = checkAlgorithm(algo);

    if(algorithmNumber == -1){
        MITK_ERROR << "Incorrect algorithm identifier. See help below: ";
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }


    // Default values for optional arguments
    std::string inImage = "";
    std::string inPts = "";
    std::string inElem = "";
    std::string inBp = "";
    std::string outFilename = "";
    auto verbose = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("image")){
        inImage = us::any_cast<std::string>(parsedArgs["image"]);
    }
    if (parsedArgs.end() != parsedArgs.find("points")){
        inPts = us::any_cast<std::string>(parsedArgs["points"]);
    }
    if (parsedArgs.end() != parsedArgs.find("elements")){
        inElem = us::any_cast<std::string>(parsedArgs["elements"]);
    }
    if (parsedArgs.end() != parsedArgs.find("bloodpool")){
        inBp = us::any_cast<std::string>(parsedArgs["bloodpool"]);
    }
    if (parsedArgs.end() != parsedArgs.find("output")){
        outFilename = us::any_cast<std::string>(parsedArgs["output"]);
    }
    if (parsedArgs.end() != parsedArgs.find("verbose")){
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }

    if(!checkAlgorithmParameters(algorithmNumber, inImage, inPts, inElem, inBp)){
        MITK_ERROR << "Exiting due to incorrect algorithm parameters.";
        return EXIT_FAILURE;
    }

    try{
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";
        QString imagePath = QString::fromStdString(inImage);
        QString pointPath = QString::fromStdString(inPts);
        QString elemPath = QString::fromStdString(inElem);
        QString bpPath = QString::fromStdString(inBp);
        QString outputPath = QString::fromStdString(outFilename);
        bool success;
        if(algorithmNumber==1){ // original_coords
            MITK_INFO << "ORIGINAL COORDINATES";
            success = originalCoordinates(imagePath, pointPath, outputPath, verbose);
        } else if(algorithmNumber==2){ // calc_cog
            MITK_INFO << "CALCULATE CENTRE OF GRAVITY";
            success = calculateCentreOfGravity(pointPath, elemPath, outputPath, verbose);
        } else if(algorithmNumber==3){ // region_map
            MITK_INFO << "REGION MAPPING";
            success = regionMapping(bpPath, pointPath, elemPath, outputPath, verbose);
        } else if(algorithmNumber==4){
            QFileInfo fi(imagePath);
            QFileInfo fm(pointPath);
            QString mybasename = fm.baseName();

            QString path2files = fi.absolutePath() + mitk::IOUtil::GetDirectorySeparator();
            QString shiftPts = path2files + mybasename + "_shift.pts";

            if(originalCoordinates(imagePath, pointPath, shiftPts, verbose)){
                MITK_INFO(verbose) << "[Step1] Shift to image coordinates... done.";

                QString cogPts = path2files + mybasename + "_COG.pts";
                if(calculateCentreOfGravity(shiftPts, elemPath, cogPts, verbose)){
                    MITK_INFO(verbose) << "[Step2] Calculation of COG ... done.";

                    success = regionMapping(bpPath, cogPts, elemPath, outputPath, verbose);

                    MITK_INFO(verbose && success) << "[Step3] Region mapping ... done.";
                    MITK_ERROR(!success) << "Error calculating Region Mapping.";
                } else {
                    MITK_ERROR << "Error calculating Centre of Gravity.";
                    success = false;
                }
            } else{
                MITK_ERROR << "Error calculating original coordinates.";
                success = false;
            }
        }
        MITK_INFO(success) << "Successful program";
        MITK_INFO(verbose) << "Goodbye!";

    }
    catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    }
    catch(...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }

}

bool originalCoordinates(QString imagePath, QString pointPath, QString outputPath, bool v){
    bool success = false;
    if(QFileInfo::exists(imagePath) && QFileInfo::exists(pointPath)){
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(imagePath.toStdString());
        ImageType::Pointer itkInput = ImageType::New();
        ImageType::PointType origin;
        mitk::CastToItkImage(image, itkInput);
        origin = image->GetOrigin();
        bool outputToConsole = outputPath.isEmpty();

        // MITK_INFO(v) << "Reorienting image to RAI";
        // typedef itk::OrientImageFilter<ImageType,ImageType> OrientImageFilterType;
        // OrientImageFilterType::Pointer orienter = OrientImageFilterType::New();
        // orienter->UseImageDirectionOn();
        // orienter->SetDesiredCoordinateOrientationToAxial(); // RAI
        // orienter->SetInput(itkInput);
        // orienter->Update();
        //
        // ImageType::SpacingType spacing = orienter->GetOutput()->GetSpacing();
        // ImageType::DirectionType transMat = orienter->GetOutput()->GetDirection();
        // ImageType::PointType origin = orienter->GetOutput()->GetOrigin();

        // std::cout << "ORIGIN: " << origin << '\n';
        // std::cout << "DIRECTION: " << transMat << '\n';
        // std::cerr << "SPACING: " << spacing << '\n';

        // origin[0] = -1.0*origin[0];
        // origin[1] = -1.0*origin[1];
        //
        // transMat[0][0] = -1.0*transMat[0][0];
        // transMat[1][1] = -1.0*transMat[1][1];

        std::ifstream pointFileRead;

        int nPts;
        pointFileRead.open(pointPath.toStdString());
        pointFileRead >> nPts;
        
        double x,y,z;
        double scaling = 1000;

        std::ofstream outputFileWrite;
        if(!outputToConsole){
            outputFileWrite.open(outputPath.toStdString());
            outputFileWrite << nPts << std::endl;
        } else{
            std::cout << nPts << std::endl;
        }

        // std::vector<double> pointsVector(nPts*3, 0.0);
        // std::vector<double> onePt(3,0.0);

        double xt=0.0, yt=0.0, zt=0.0;
        for(int iPt=0; iPt < nPts; iPt++){
            pointFileRead >> x;
            pointFileRead >> y;
            pointFileRead >> z;
            if (pointFileRead.eof()) {
                cerr << " WARNING!: File ended prematurely " << endl;
                break;
            }

            xt = x+(origin[0]*scaling);
            yt = y+(origin[1]*scaling);
            zt = z+(origin[2]*scaling);

            // for(int ix=0; ix<3; ix++){
                // pointFileRead >> pointsVector[iPt + ix*nPts];
                // pointFileRead >> onePt[ix];
                // onePt[ix] /= 1000;
            // }
            // xt = transMat[0][0] * onePt[0] + transMat[0][1] * onePt[1] + transMat[0][2] * onePt[2] + origin[0];
            // yt = transMat[1][0] * onePt[0] + transMat[1][1] * onePt[1] + transMat[1][2] * onePt[2] + origin[1];
            // zt = transMat[2][0] * onePt[0] + transMat[2][1] * onePt[1] + transMat[2][2] * onePt[2] + origin[2];
            //
            // xt *= 1000;
            // yt *= 1000;
            // zt *= 1000;


            // onePt.clear();

            if(!outputToConsole){
                outputFileWrite << std::fixed << xt << " ";
                outputFileWrite << std::fixed << yt << " ";
                outputFileWrite << std::fixed << zt << std::endl;
            } else{
                std::cout << std::fixed << xt << std::fixed << yt << std::fixed << zt << std::endl;
            }
        }
        pointFileRead.close();

        success = true;

    } else{
        MITK_ERROR(QFileInfo::exists(imagePath)) << ("Could not read file" + imagePath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("Could not read file" + pointPath).toStdString();
    }
    return success;
}

bool calculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath, bool v){
    bool success = false;
    if(QFileInfo::exists(elemPath) && QFileInfo::exists(pointPath)){
        std::ifstream pointFileRead;
        std::ifstream elemFileRead;
        int nPts, nElem;
        bool outputToConsole = outputPath.isEmpty();

        pointFileRead.open(pointPath.toStdString());
        pointFileRead >> nPts;

        std::vector<double> pointsVector(nPts*3, 0.0);
        for(int iPt=0; iPt < nPts; iPt++){
            for(int ix=0; ix<3; ix++){
                pointFileRead >> pointsVector[iPt + ix*nPts];
            }
        }
        pointFileRead.close();

        elemFileRead.open(elemPath.toStdString());
        elemFileRead >> nElem;

        std::ofstream outputFileWrite;
        if(!outputToConsole){
            outputFileWrite.open(outputPath.toStdString());
            outputFileWrite << nElem << " 3"<< std::endl;
        } else{
            std::cout << nElem <<" 3"<< std::endl;
        }

        std::string type;
        for(int iElem=0; iElem < nElem; iElem++){
            int p[4], region;
            double x=0.0, y=0.0, z=0.0;

            elemFileRead >> type;
            elemFileRead >> p[0];
            elemFileRead >> p[1];
            elemFileRead >> p[2];
            elemFileRead >> p[3];
            elemFileRead >> region;

            for(int ix=0; ix<4; ix++){
                x += pointsVector.at(p[ix]+0);
                y += pointsVector.at(p[ix]+1);
                z += pointsVector.at(p[ix]+2);
            }

            x /= 4.0 *1000;
            y /= 4.0 *1000;
            z /= 4.0 *1000;

            if(!outputToConsole){
                outputFileWrite << std::fixed << x << std::endl;
                outputFileWrite << std::fixed << y << std::endl;
                outputFileWrite << std::fixed << z << std::endl;
            } else{
                std::cout << std::fixed << x << std::endl;
                std::cout << std::fixed << y << std::endl;
                std::cout << std::fixed << z << std::endl;
            }
        }

        if(!outputToConsole){
            outputFileWrite.close();
        }
        success = true;
        MITK_INFO(v) << "Completed input .elem file";


    } else{
        MITK_ERROR(QFileInfo::exists(elemPath)) << ("Could not read file" + elemPath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("Could not read file" + pointPath).toStdString();
    }

    return success;

}

bool regionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath, bool v){
    bool success = false;
    if(QFileInfo::exists(bpPath) && QFileInfo::exists(pointPath) && QFileInfo::exists(elemPath)){
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(bpPath.toStdString());
        ImageType::Pointer itkInput = ImageType::New();
        mitk::CastToItkImage(image, itkInput);
        bool outputToConsole = outputPath.isEmpty();

        ImageType::SpacingType spacing = itkInput->GetSpacing();
        ImageType::RegionType imageRegion = itkInput->GetLargestPossibleRegion();
        ImageType::SizeType size = imageRegion.GetSize();
        int minIx, minJx, minKx, maxIx, maxJx, maxKx;
        minIx=0; minJx=0; minKx=0;
        maxIx=size[0]-1;
        maxJx=size[1]-1;
        maxKx=size[2]-1;

        std::ofstream outputFileWrite;

        std::ifstream cogFileRead, elemFileRead;
        cogFileRead.open(pointPath.toStdString());

        int nElemCOG, dim, count;
        double x, y, z;

        cogFileRead >> nElemCOG;
        cogFileRead >> dim;

        if(!outputToConsole){
            outputFileWrite.open(outputPath.toStdString());
            outputFileWrite << nElemCOG << std::endl;
        } else{
            std::cout << nElemCOG << std::endl;
        }

        MITK_INFO(v) << ("Number of elements (COG file):" + QString::number(nElemCOG)).toStdString();
        MITK_INFO(v) << ("Dimension= " + QString::number(dim)).toStdString();

        elemFileRead.open(elemPath.toStdString());
        int nElem;

        elemFileRead >> nElem;
        if(nElem != nElemCOG){
            MITK_ERROR << "Number of elements in files are not consistent.";
        }

        count = 0;
        std::string type;
        int p[4], region, newRegion;
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
            elemFileRead >> p[0];
            elemFileRead >> p[1];
            elemFileRead >> p[2];
            elemFileRead >> p[3];
            elemFileRead >> region;

            newRegion = inRegion(itkInput, x, y, z);
            if(newRegion > 0 ){
                region = newRegion;
                newRegionCount++;
            }

            if(!outputToConsole){
                outputFileWrite << type << " ";
                outputFileWrite << p[0] << " ";
                outputFileWrite << p[1] << " ";
                outputFileWrite << p[2] << " ";
                outputFileWrite << p[3] << " ";
                outputFileWrite << region << std::endl;
            } else{
                std::cout << type << " " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << region << std::endl;
            }
            count++;
        }

        if(!outputToConsole){
            outputFileWrite.close();
        }
        success = true;

        MITK_INFO << ("Number of element COG read: " + QString::number(count)).toStdString();
        MITK_INFO << ("Number of new regions determined: " + QString::number(newRegionCount)).toStdString();

        if ( minIx < 0 || minJx < 0 || minKx < 0 ) {
            MITK_WARN << "WARNING: The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region.";
            MITK_WARN << "If scar does lie in this region, then you need to pad the image at the start by (in pixel space):";
            MITK_WARN << (QString::number(-minIx) + " " + QString::number(-minJx) + " " + QString::number(-minKx)).toStdString();
            MITK_WARN << "And add the following transformation to the TransMatFile (in geometric space):";
            MITK_WARN << (QString::number(-minIx*spacing[0]) + " " + QString::number(-minJx*spacing[1]) + " " + QString::number(-minKx*spacing[2])).toStdString();
            success=false;
        }
        if ( maxIx > int(size[0]-1) || maxJx > int(size[1]-1) || maxKx > int(size[2]-1) ) {
            MITK_WARN << "WARNING: The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region.";
            MITK_WARN << "If scar does lie in this region, then you need to pad the image at the end by (in pixel space):";
            MITK_WARN << (QString::number(maxIx-(size[0]-1)) + " " + QString::number(maxJx-(size[1]-1)) + " " + QString::number(maxKx-(size[2]-1))).toStdString();
            MITK_WARN << "No need to change TransMatFile";
            success=false;
        }

    } else{
        MITK_ERROR(QFileInfo::exists(bpPath)) << ("File does not exist: " + bpPath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("File does not exist: " + pointPath).toStdString();
    }

    return success;
}

// helper functions
int inRegion(ImageType::Pointer image, double x, double y, double z){
    ImageType::SpacingType spacing = image->GetSpacing();
    ImageType::RegionType region = image->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    int minIx, minJx, minKx, maxIx, maxJx, maxKx;
    minIx=0; minJx=0; minKx=0;
    maxIx=size[0]-1;
    maxJx=size[1]-1;
    maxKx=size[2]-1;

    int ix, jx, kx;
    ix = static_cast<int>(x/spacing[0]);
    jx = static_cast<int>(y/spacing[1]);
    kx = static_cast<int>(z/spacing[2]);

    if(ix<0 || jx<0 || kx<0){
        minIx = std::min(ix, minIx);
        minJx = std::min(jx, minJx);
        minKx = std::min(kx, minKx);
        return 0;
    } else if(ix > maxIx || jx > maxJx || kx > maxKx){
        maxIx = std::max(ix, maxIx);
        maxJx = std::max(jx, maxJx);
        maxKx = std::max(kx, maxKx);
        return 0;
    } else{
        IndexType index={{ix, jx, kx}};
        return image->GetPixel(index);
    }
}

int checkAlgorithm(std::string algo){
    QString algorithm = QString::fromStdString(algo);

    MITK_INFO << "Algorithm input: " + algo;

    int algorithmNumber = -1;

    if(algo.compare("original_coords")==0 || algo.compare("coords")==0){
        algorithmNumber = 1;
    }
    if(algo.compare("calc_cog")==0 || algo.compare("cog")==0){
        algorithmNumber = 2;
    }
    if(algo.compare("region_map")==0 || algo.compare("region")==0){
        algorithmNumber = 3;
    }
    if(algo.compare("all")==0){
        algorithmNumber = 4;
    }

    return algorithmNumber;
}

bool checkAlgorithmParameters(int algorithmNumber, std::string inImage, std::string inPts, std::string inElem, std::string inBp){

    bool parametersCorrect = false;
    bool imageEmpty = inImage.empty();
    bool ptsEmpty = inPts.empty();
    bool elemEmpty = inElem.empty();
    bool bpEmpty = inBp.empty();

    if(algorithmNumber == 1 && (imageEmpty || ptsEmpty || bpEmpty)){
        MITK_ERROR << "METHOD: original_coords - Incorrect parameters:";
        MITK_INFO << "\nUsage:\n\t./MitkCemrgCarpUtils --method original_coords -i /path/to/image.nii -pts /path/to/points.pts";
    } else if(algorithmNumber == 2 && (ptsEmpty || elemEmpty)){
        MITK_ERROR << "METHOD: calc_cog - Incorrect parameters:";
        MITK_INFO << "\nUsage:\n\t./MitkCemrgCarpUtils --method calc_cog -pts /path/to/points.pts -elem /path/to/elements.elem";
    } else if(algorithmNumber == 3 && (imageEmpty || ptsEmpty || elemEmpty)){
        MITK_ERROR << "METHOD: region_map - Incorrect parameters:";
        MITK_INFO << "\nUsage:\n\t./MitkCemrgCarpUtils --method region_map -pts /path/to/points.pts -elem /path/to/elements.elem -bp /path/to/bloodpool.nii";
    } if(algorithmNumber == 4 && (imageEmpty || ptsEmpty || bpEmpty)){
        MITK_ERROR << "ALL METHODS workflow: Incorrect parameters";
        MITK_INFO << "\nUsage:\n\t./MitkCemrgCarpUtils --method all -i /path/to/image.nii -pts /path/to/points.pts -elem /path/to/elements.elem -bp /path/to/bloodpool.nii";
    } else{
        parametersCorrect = true;
    }

    MITK_INFO(imageEmpty) << "Specify the image path with flag --image (or -i).";
    MITK_INFO(ptsEmpty) << "Specify the points file path with flag --points (or -pts).";
    MITK_INFO(elemEmpty) << "Specify the elements file path with flag --elements (or -elem).";
    MITK_INFO(bpEmpty) << "Specify the bloodpool parameter with flag --bloodpool (or -bp).";

    return parametersCorrect;
}
