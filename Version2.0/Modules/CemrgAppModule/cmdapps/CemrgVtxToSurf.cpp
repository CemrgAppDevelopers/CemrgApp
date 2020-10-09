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
#include <QDir>
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

void saveVectorToFile(QString outputfilename, std::vector<int> & indicesList);
// std::set<int> & vector2set(std::vector<int> v);

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("EASI processing");
    parser.setTitle("Extract Surface from Labels");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Extract Surface from Labels CemrgApp Command-line App.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::InputFile,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
                "mesh", "msh", mitkCommandLineParser::String,
                "Mesh path.", "Path to mesh (.elem) file.");
    parser.addArgument(
                "lv_labels", "lv", mitkCommandLineParser::String,
                "Labels corresponding to LV: epi,endo,base", "Syntax: -lv 1,10,2");
    parser.addArgument(
                "rv_labels", "rv", mitkCommandLineParser::String,
                "Labels corresponding to RV: epi,endo,base", "Syntax: -rv 5,30,4");
    parser.addArgument(
                "apex_labels", "apex", mitkCommandLineParser::String,
                "Label corresponding to apex", "Syntax: -apex 3");
    parser.addArgument( // optional
                "verbose", "v", mitkCommandLineParser::Bool,
                "Verbose Output", "Whether to produce verbose output");
    parser.addArgument( // optional
                "cavity-mesh", "cav", mitkCommandLineParser::Bool,
                "Use <meshName>_Cav.elem", "Change <meshName>.elem when reading the mesh.");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty()){
        return EXIT_FAILURE;
    }

    if (parsedArgs["mesh"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    auto meshfile = us::any_cast<std::string>(parsedArgs["mesh"]);

    // Default values for optional arguments
    std::string labs_LV = "1,10,2";
    std::string labs_RV = "5,30,4";
    std::string labs_Apex = "3";
    auto cavFlag = false;
    auto verbose = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("lv_labels")){
        labs_LV = us::any_cast<std::string>(parsedArgs["lv_labels"]);
    }
    if (parsedArgs.end() != parsedArgs.find("rv_labels")){
        labs_RV = us::any_cast<std::string>(parsedArgs["rv_labels"]);
    }
    if (parsedArgs.end() != parsedArgs.find("apex_labels")){
        labs_Apex = us::any_cast<std::string>(parsedArgs["apex_labels"]);
    }
    if (parsedArgs.end() != parsedArgs.find("cavity-mesh")){
        cavFlag = us::any_cast<bool>(parsedArgs["cavity-mesh"]);
    }
    if (parsedArgs.end() != parsedArgs.find("verbose")){
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }

    try{
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";
        QFileInfo fi(QString::fromStdString(meshfile));

        MITK_INFO(verbose) << "[Step1] Extracting labels";
        QString labelsLeft = QString::fromStdString(labs_LV);
        QString labelsRight = QString::fromStdString(labs_RV);
        QString labelsApex = QString::fromStdString(labs_Apex);

        QStringList labelsLeftList = labelsLeft.split(QLatin1Char(','));
        QStringList labelsRightList = labelsRight.split(QLatin1Char(','));

        int lvepi = labelsLeftList.at(0).toInt();
        int lvendo = labelsLeftList.at(1).toInt();
        int lvbase = labelsLeftList.at(2).toInt();

        int rvepi = labelsRightList.at(0).toInt();
        int rvendo = labelsRightList.at(1).toInt();
        int rvbase = labelsRightList.at(2).toInt();

        int apex = labelsApex.toInt();

        std::vector<int> lvepiElemList, lvendoElemList, lvbaseElemList;
        std::vector<int> rvepiElemList, rvendoElemList, rvbaseElemList, apexElemList;

        std::vector<int> labelsToFile;
        labelsToFile.push_back(lvendo);
        labelsToFile.push_back(lvepi);
        labelsToFile.push_back(lvbase);
        labelsToFile.push_back(rvendo);
        labelsToFile.push_back(rvepi);
        labelsToFile.push_back(rvbase);
        labelsToFile.push_back(apex);

        QString path2files = fi.absolutePath() + mitk::IOUtil::GetDirectorySeparator();
        QString basename = fi.baseName();
        QString readInMeshFile = cavFlag ? basename + "_Cav.elem" : fi.completeBaseName();
        QString readInPtsFile = basename + ".pts";
        QString outputBase = path2files + basename + "_";

        QStringList outputFiles;
        outputFiles << (outputBase+"LVendo.surf.vtx");
        outputFiles << (outputBase+"LVepi.surf.vtx");
        outputFiles << (outputBase+"LVbase.surf.vtx");
        outputFiles << (outputBase+"RVendo.surf.vtx");
        outputFiles << (outputBase+"RVepi.surf.vtx");
        outputFiles << (outputBase+"RVbase.surf.vtx");
        outputFiles << (outputBase+"apex.surf.vtx");

        std::ifstream elemFileRead, ptsFileRead;
        std::ofstream outputFileWrite;

        int nElem, nPts;

        MITK_INFO(verbose) << ("[...] Reading in mesh file: " + path2files + readInMeshFile).toStdString();
        ptsFileRead.open((path2files + readInPtsFile).toStdString());
        elemFileRead.open((path2files + readInMeshFile).toStdString());
        ptsFileRead >> nPts;
        elemFileRead >> nElem;
        ptsFileRead.close();

        int totalTags = 7;
        std::vector<int> ptsVector(nPts*totalTags, 0);

        MITK_INFO(verbose) << ("[...] Number of elements in mesh: " + QString::number(nElem)).toStdString();
        MITK_INFO(verbose) << ("[...] Number of points in mesh: " + QString::number(nPts)).toStdString();

        MITK_INFO(verbose) << "[Step2] Read in regions from elements file.";
        std::string type;
        for(int iElem=0; iElem < nElem; iElem++){
            int p[4], currentTag;
            elemFileRead >> type;
            elemFileRead >> p[0];
            elemFileRead >> p[1];
            elemFileRead >> p[2];
            elemFileRead >> p[3];
            elemFileRead >> currentTag;

            for(int ix=0; ix<4; ix++){
                int it=0;
                try{
                    while(it<7 && ptsVector.at(it + totalTags*p[ix]) != 0){
                        it++;
                    }
                    if(it<7){
                        ptsVector.at(it + totalTags*p[ix]) = currentTag;
                    }
                } catch(...){
                    MITK_WARN << ("Size: " + QString::number(ptsVector.size())).toStdString();
                    MITK_WARN << ("it: " + QString::number(it) + " p[ix]:" + QString::number(p[ix])).toStdString();
                }
            }

            if(verbose && iElem%500000==0){
                MITK_INFO << ("iElem: " + QString::number(iElem)).toStdString();
            }
        }
        MITK_INFO(verbose) << "[...] Finished reading regions";

        elemFileRead.close();

        QDir dir(path2files);
        QStringList filters;
        filters << "*.vtx";
        dir.setNameFilters(filters);
        QStringList surfFilesList = dir.entryList();

        MITK_INFO(verbose) << "[Step3] Selecting indices based on elements' regions";
        for(int ix=0; ix< surfFilesList.size(); ix++){
            // MITK_INFO << (basename + " " + surfFilesList.at(ix)).toStdString();
            if(surfFilesList.at(ix).contains(basename, Qt::CaseSensitive)){
                std::ifstream surfFile;
                QString path2surffile = path2files + surfFilesList.at(ix);
                surfFile.open(path2surffile.toStdString());

                MITK_INFO << path2surffile.toStdString();

                int nElemInSurf;
                surfFile >> nElemInSurf;
                MITK_INFO(verbose) << ("Num Points in Surface: " + QString::number(nElemInSurf)).toStdString();
                for(int jx=0; jx<nElemInSurf; jx++){
                    // for(int jx=0; jx<4; jx++){ // for testing
                    int currentIndex, currentLabel;
                    surfFile >> currentIndex;
                    for(int kx=0; kx<totalTags; kx++){
                        currentLabel = ptsVector.at(kx + currentIndex*totalTags);
                        if(jx<2){
                            std::cout << "Curr Label:" << currentLabel << '\n';
                        }

                        if(currentLabel==lvepi){
                            lvepiElemList.push_back(currentIndex);
                        } else if(currentLabel==lvbase){
                            lvbaseElemList.push_back(currentIndex);
                        } else  if(currentLabel==apex){
                            apexElemList.push_back(currentIndex);
                        } else if(currentLabel==rvbase){
                            rvbaseElemList.push_back(currentIndex);
                        } else if(currentLabel==rvepi){
                            rvepiElemList.push_back(currentIndex);
                        } else if(currentLabel==lvendo){
                            lvendoElemList.push_back(currentIndex);
                        } else if(currentLabel==rvendo){
                            rvendoElemList.push_back(currentIndex);
                        }
                    }
                }
                surfFile.close();
            }
        }



        MITK_INFO(verbose) << "[Step4] Saving to corresponding .surf.vtx files";
        for(int ix=0; ix<outputFiles.size();ix++){
            int currentLabel = labelsToFile.at(ix);
            if(currentLabel==lvepi){
                saveVectorToFile(outputFiles.at(ix), lvepiElemList);
            } else if(currentLabel==lvbase){
                saveVectorToFile(outputFiles.at(ix), lvbaseElemList);
            } else  if(currentLabel==apex){
                saveVectorToFile(outputFiles.at(ix), apexElemList);
            } else if(currentLabel==rvbase){
                saveVectorToFile(outputFiles.at(ix), rvbaseElemList);
            } else if(currentLabel==rvepi){
                saveVectorToFile(outputFiles.at(ix), rvepiElemList);
            } else if(currentLabel==lvendo){
                saveVectorToFile(outputFiles.at(ix), lvendoElemList);
            } else if(currentLabel==rvendo){
                saveVectorToFile(outputFiles.at(ix), rvendoElemList);
            }
        }

        bool success = true;

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

void saveVectorToFile(QString outputfilename, std::vector<int> & indicesList){
    MITK_INFO << ("[...] Saving file" + outputfilename).toStdString();
    std::ofstream outfile;

    outfile.open(outputfilename.toStdString());
    if(!indicesList.empty()){
        std::vector<int>::iterator ip;
        std::sort(indicesList.begin(), indicesList.end());
        ip = std::unique(indicesList.begin(), indicesList.end());
        indicesList.resize(std::distance(indicesList.begin(), ip));

        outfile << indicesList.size() << std::endl;
        outfile << "extra" << std::endl;
        for(int elemIdx : indicesList){
            outfile << elemIdx << std::endl;
        }

    } else{
        outfile << "0";
    }

    outfile.close();
}
//
// std::set<int> & vector2set(std::vector<int> v){
//     std::set<int> s(v.begin(), v.end());
//     return s;
// }
