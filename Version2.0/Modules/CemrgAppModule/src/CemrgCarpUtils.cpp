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
#include <vtkImageMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor2D.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkCenterOfMass.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkExtractVOI.h>
#include <vtkPlane.h>
#include <vtkProperty.h>
#include <vtkCutter.h>
#include <vtkCamera.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkColorTransferFunction.h>

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
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFileDialog>
#include <QTextStream>

#include "CemrgCarpUtils.h"

// basics
CemrgCarpUtils::CemrgCarpUtils(){
    elemPath = "";
    ptsPath = "" ;
    dir = "";
    bname = "";

    tol = 1e-7;

    alpha_epi = -50;
    alpha_endo = 40;
    beta_endo = -65;
    beta_epi = 25;
}

void CemrgCarpUtils::UpdateFileNames(){
    QFileInfo elfi(elemPath);
    bname = elfi.baseName();
    dir = elfi.absolutePath();
    ptsPath = dir + mitk::IOUtil::GetDirectorySeparator() + bname + ".pts";
}

void CemrgCarpUtils::ReadElementFile(int nIter){
    if(elemPath.isEmpty()){
        MITK_ERROR << "Element path name not set.";
        return;
    } else{
        std::ifstream fi(elemPath.toStdString());
        if (fi){
            fi >> nElem;
            if(nIter>0){
                MITK_INFO << ("Chaning total iterations to: " + QString::number(nIter)).toStdString();
                nElem = nIter;
            }

            el.resize(nElem*4);
            for (int ix = 0; ix < nElem; ix++) {
                std::string type;
                int region;
                fi >> type;
                for (int jx = 0; jx < 4; jx++) {
                    fi >> el[4*ix + jx];
                }
                fi >> region;
            }
            fi.close();
        } else{
            MITK_ERROR << ("Error opening file " + elemPath).toStdString();
            return;
        }
    }
}

void CemrgCarpUtils::ReadPointsFile(){
    if(ptsPath.isEmpty()){
        MITK_ERROR << "Points path name not set.";
        return;
    } else{
        std::ifstream fi(ptsPath.toStdString());
        if (fi){
            fi >> nPts;
            pts.resize(nPts*3);
            for (int ix = 0; ix < nPts; ix++) {
                for (int jx = 0; jx < 3; jx++) {
                    fi >> pts[3*ix + jx];
                }
            }
            fi.close();
        } else{
            MITK_ERROR << ("Error opening file " + ptsPath).toStdString();
            return;
        }
    }
}

void CemrgCarpUtils::ReadLaplaceSolvesFile(QString dataFilePath, int whichFile){
    if(whichFile==0){ // apex
        apexPts = ReadInFile(dataFilePath);
    } else if(whichFile==1){ //epi
        epiPts = ReadInFile(dataFilePath);
    } else if(whichFile==2){ // lv
        lvPts = ReadInFile(dataFilePath);
    } else if(whichFile==3){ // rv
        rvPts = ReadInFile(dataFilePath);
    }
}


void CemrgCarpUtils::ReadElementGradientFile(QString gradientPath, int whichFile){
    if(whichFile==0){ // apex
        apexGrad = ReadInGradient(gradientPath);
    } else if(whichFile==1){ //epi
        epiGrad = ReadInGradient(gradientPath);
    } else if(whichFile==2){ // lv
        lvGrad = ReadInGradient(gradientPath);
    } else if(whichFile==3){ // rv
        rvGrad = ReadInGradient(gradientPath);
    }
}

double CemrgCarpUtils::DataCentre(std::vector<double> fullDatVector, std::vector<int> elemIdx){
    double result = 0.0;
    int N = elemIdx.size();
    if (N>0){
        for (int iv = 0; iv < N; iv++) {
            result += fullDatVector[elemIdx[iv]];
        }
        result /= N;
    }
    return result;
}

std::vector<double> CemrgCarpUtils::ReadInFile(QString dataFilePath){
    std::ifstream fi(dataFilePath.toStdString());
    std::vector<double> result;
    if(fi){
        result.resize(nPts);
        for (int ix = 0; ix < nPts; ix++) {
            fi >> result[ix];
        }
        fi.close();
    } else{
        MITK_ERROR << ("Error opening file " + dataFilePath).toStdString();
    }
    return result;
}

std::vector<double> CemrgCarpUtils::ReadInGradient(QString gradientPath){
    std::ifstream fi(gradientPath.toStdString());
    std::vector<double> result;
    if(fi){
        result.resize(nElem*3);
        for (int ix = 0; ix < nElem; ix++) {
            for (int jx = 0; jx < 3; jx++) {
                fi >> result[3*ix + jx];
            }
        }
        fi.close();
    } else{
        MITK_ERROR << ("Error opening file " + gradientPath).toStdString();
    }
    return result;
}

// Fibres calculation functions
Quaternion CemrgCarpUtils::QuaternionMult(const Quaternion a, const Quaternion b){
    Quaternion result;
    result.s = a.s*b.s - a.i*b.i - a.j*b.j - a.k*b.k;
    result.i = a.s*b.i + a.i*b.s + a.j*b.k - a.k*b.j;
    result.j = a.s*b.j - a.i*b.k + a.j*b.s + a.k*b.i;
    result.k = a.s*b.k + a.i*b.j - a.j*b.i + a.k*b.s;
    return result;
}

Quaternion CemrgCarpUtils::QuaternionNormalise(const Quaternion q){
    double norm;
    norm = sqrt(QuaternionDot(q, q));
    Quaternion result;
    result.s = q.s / norm;
    result.i = q.i / norm;
    result.j = q.j / norm;
    result.k = q.k / norm;
    return result;
}

Quaternion CemrgCarpUtils::RotationToQuaternion(double M[3][3]){
    Quaternion q;
    double miniTol = 1e-12;

    double w2 = (1+M[0][0]+M[1][1]+M[2][2])/4;
    if (w2 > miniTol) {
      q.s = sqrt(w2);
      q.i = (M[1][2]-M[2][1])/4/q.s;
      q.j = (M[2][0]-M[0][2])/4/q.s;
      q.k = (M[0][1]-M[1][0])/4/q.s;
    } else {
      q.s = 0;
      double x2 = -(M[1][1]+M[2][2])/2;
      if (x2 > miniTol) {
        q.i = sqrt(x2);
        q.j = M[0][1]/2/q.i;
        q.k = M[0][2]/2/q.i;
      } else {
        q.i = 0;
        double y2 = (1-M[2][2])/2;
        if (y2 > miniTol) {
          q.j = sqrt(y2);
          q.k = M[1][2]/2/q.j;
        } else {
          q.j = 0;
          q.k = 1;
        }
      }
    }
    return q;
}

void CemrgCarpUtils::QuaternionToRotation(const Quaternion q, double M[3][3]){
    double w=q.s;
    double x=q.i;
    double y=q.j;
    double z=q.k;
    M[0][0] = w*w+x*x-y*y-z*z;
    M[0][1] = 2*x*y+2*w*z;
    M[0][2] = 2*x*z-2*w*y;

    M[1][0] = 2*x*y-2*w*z;
    M[1][1] = w*w-x*x+y*y-z*z;
    M[1][2] = 2*y*z+2*w*x;

    M[2][0] = 2*x*z+2*w*y;
    M[2][1] = 2*y*z-2*w*x;
    M[2][2] = w*w-x*x-y*y+z*z;
}

double CemrgCarpUtils::DotProd(const double a[3], const double b[3]){
    double sum=0;
    for (int ii=0; ii<3; ii++) {
      sum += a[ii]*b[ii];
    }
    return sum;
}

void CemrgCarpUtils::NormaliseVector(double v[3]){
    double mag = sqrt(DotProd(v,v));
    for (int ii=0; ii<3; ii++) {
      v[ii] /= mag;
    }
}

void CemrgCarpUtils::CrossProd(const double a[3], const double b[3], double c[3]){
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = b[0]*a[2] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

Quaternion CemrgCarpUtils::AxisSystem(const std::vector<double> apex, std::vector<double> surface_inward){
    //Compute the reference frame
    double M[3][3];

    double e0[3];
    for (int ii=0; ii<3; ii++) {
      e0[ii] = apex[ii];
    }
    NormaliseVector(e0);
    for (int ii=0; ii<3; ii++) {
      M[ii][0] = e0[ii];
    }
    double e2[3];
    double surf_array[3] = {surface_inward[0], surface_inward[1], surface_inward[2]};
    double dotp = DotProd(e0, surf_array);
    for (int ii=0; ii<3; ii++) {
      e2[ii] = surface_inward[ii] - dotp*e0[ii];
    }
    NormaliseVector(e2);
    for (int ii=0; ii<3; ii++) {
      M[ii][2] = e2[ii];
    }
    double e1[3];
    CrossProd(e2,e0,e1);
    NormaliseVector(e1);
    for (int ii=0; ii<3; ii++) {
      M[ii][1] = e1[ii];
    }

    return RotationToQuaternion(M);
}

Quaternion CemrgCarpUtils::RotateAxis(const Quaternion q, const double alpha, const double beta){
    Quaternion a = {cos(-alpha/2),0,            0, sin(-alpha/2)};
    Quaternion b = {cos(-beta/2) ,sin(-beta/2), 0, 0            };
    return QuaternionMult(QuaternionMult(b,a),q);
}
Quaternion CemrgCarpUtils::Slerp(const double a_coeff, const Quaternion a, const double b_coeff, const Quaternion b){
    Quaternion R[4];
    //R0 = a
    R[0] = a;
    //R1 = i*a
    R[1].s = -a.i;
    R[1].i =  a.s;
    R[1].j = -a.k;
    R[1].k =  a.j;
    //R2 = j*a
    R[2].s = -a.j;
    R[2].i =  a.k;
    R[2].j =  a.s;
    R[2].k = -a.i;
    //R3 = k*a
    R[3].s = -a.k;
    R[3].i = -a.j;
    R[3].j =  a.i;
    R[3].k =  a.s;

    //find the quaternion from a canidates that is closest to b.
    double max_dot_p=0;
    int max_quat_index=-1;
    for (int ii=0; ii<4; ii++) {
      double dot_p = QuaternionDot(R[ii],b);
      if (dot_p < 0) {
        dot_p = -dot_p;
        R[ii].s *= -1;
        R[ii].i *= -1;
        R[ii].j *= -1;
        R[ii].k *= -1;
      }
      if (dot_p > max_dot_p) {
        max_quat_index=ii;
        max_dot_p = dot_p;
      }
    }
    Quaternion max_quat = R[max_quat_index];

    Quaternion result;
    if (max_dot_p > 1-tol) {
      //in this case, theta == sin_theta == 0.  Quaterions are so close, so just pick one.
      result = b;
    } else {
      //compute the interpolation.
      double theta = acos(max_dot_p);
      double sin_theta = sqrt(1-max_dot_p*max_dot_p);

      double adjusted_a_coeff = sin(theta*a_coeff)/sin_theta;
      double adjusted_b_coeff = sin(theta*b_coeff)/sin_theta;
      result.s = max_quat.s*adjusted_a_coeff + b.s*adjusted_b_coeff;
      result.i = max_quat.i*adjusted_a_coeff + b.i*adjusted_b_coeff;
      result.j = max_quat.j*adjusted_a_coeff + b.j*adjusted_b_coeff;
      result.k = max_quat.k*adjusted_a_coeff + b.k*adjusted_b_coeff;
    }
    return result;
}

Quaternion CemrgCarpUtils::BiVInterpolate(Quaternion epi_axis, Quaternion lv_axis, Quaternion rv_axis,
                           double epi, double lv, double rv){
    Quaternion q;

    if (lv < tol && rv < tol) {
        // Then we are on the epicardium - just compute epicardial fibres
        q = RotateAxis(epi_axis, (alpha_epi + 90.0)*M_PI/180.0, beta_epi*M_PI/180.0);

    } else {
        // as = alpha_septum, bs=beta_septum
        double septum = lv + rv;
        double as = ((90.0 - alpha_endo) * (lv/septum) + (90.0 + alpha_endo) * (rv/septum)) * M_PI/180.0;
        double bs  = ((-1.0*beta_endo) * (lv/septum) + beta_endo * (rv/septum)) * M_PI/180.0;
        Quaternion qlv;
        if (lv >= tol) {
            qlv = RotateAxis(lv_axis, as, bs);
        }

        Quaternion qrv;
        if (rv >= tol) {
            qrv = RotateAxis(rv_axis, -1.0*as, bs);
        }

        Quaternion qseptum;
        if (lv >= tol && rv < tol) {
           qseptum = qlv;
        } else if (lv < tol && rv >= tol) {
           qseptum = qrv;
        } else {
           qseptum = Slerp(lv/septum, qlv, rv/septum, qrv);
        }
        if (epi < tol) {
           q = qseptum;
        } else {
            // aw = alpha_wall, bw=beta_wall
            double aw = ((alpha_epi  + 90.0)*epi + (alpha_endo + 90.0)*septum) * M_PI/180.0;
            double bw  = (beta_epi*epi + beta_endo*septum) * M_PI/180.0;
            Quaternion qepi = RotateAxis(epi_axis, aw, bw);
            q = Slerp(epi, qepi, septum, qseptum);
        }
    }

    return q;
}

void CemrgCarpUtils::DefineFibres(){
    QString outPath = dir + mitk::IOUtil::GetDirectorySeparator() + bname+"_fibres.lon";
    std::ofstream fo(outPath.toStdString());
    fo << "2" << std::endl;

    std::vector<double> epi_gradient, apex_gradient, lv_gradient, rv_gradient;
    for (int qx = 0; qx < nElem; qx++) {
        std::vector<int> elemIdx(4);
        for (int ix = 0; ix < 4; ix++) {
            elemIdx[ix] = el[4*qx + ix];
        }

        double lv = DataCentre(lvPts, elemIdx);
        double rv = DataCentre(rvPts, elemIdx);
        double epi = DataCentre(rvPts, elemIdx);

        for (int jx = 0; jx < 3; jx++) {
            epi_gradient[jx] = epiGrad[3*qx + jx];
            apex_gradient[jx] = apexGrad[3*qx + jx];
            lv_gradient[jx] = lvGrad[3*qx + jx];
            rv_gradient[jx] = rvGrad[3*qx + jx];
        }

        Quaternion epi_axis, lv_axis, rv_axis;
        if (epi >= tol){
          epi_axis = AxisSystem(apex_gradient, epi_gradient);
        }

        if (lv >= tol){
          lv_axis = AxisSystem(apex_gradient, lv_gradient);
        }

        if (rv >= tol){
          rv_axis = AxisSystem(apex_gradient, rv_gradient);
        }

        Quaternion q;
        q = BiVInterpolate(epi_axis, lv_axis, rv_axis, epi, lv, rv);

        double R[3][3];
        int precision = 8;
        QuaternionToRotation(q, R);
        fo << std::fixed << std::setprecision(precision) << R[0][0] << " "<< R[1][0] << " "<< R[2][0] << " "<< R[0][2] << " "<< R[1][2] << " "<< R[2][2] << std::endl;
    }

    fo.close();
}

// static functions
void CemrgCarpUtils::OriginalCoordinates(QString imagePath, QString pointPath, QString outputPath, double scaling){
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

void CemrgCarpUtils::CalculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath){
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

            outputFileWrite << std::fixed << std::setprecision(6) << x << std::endl;
            outputFileWrite << std::fixed << std::setprecision(6) << y << std::endl;
            outputFileWrite << std::fixed << std::setprecision(6) << z << std::endl;

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

void CemrgCarpUtils::RegionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath){
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
            MITK_WARN << "The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region.";
            std::cerr << "If scar does lie in this region, then you need to pad the image at the start by (in pixel space):"<< std::endl;
            std::cerr << -min_i << " " << -min_j << " " -min_k<< std::endl;
            std::cerr << "And add the following transformation to the TransMatFile (in geometric space):"<< std::endl;
            std::cerr << -min_i*spacing[0] << " " << -min_j*spacing[1] << " " << -min_k*spacing[2]<< std::endl;
        }
        if ( max_i > int(size[0]-1) || max_j > int(size[1]-1) || max_k > int(size[2]-1) ) {
            MITK_WARN << "The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region.";
            std::cerr << "If scar does lie in this region, then you need to pad the image at the end by (in pixel space):"<< std::endl;
            std::cerr << max_i-(size[0]-1) << " " << max_j-(size[1]-1) << " " << max_k-(size[2]-1)<< std::endl;
            std::cerr << "No need to change TransMatFile"<< std::endl;
        }

    } else{
        MITK_ERROR(QFileInfo::exists(bpPath)) << ("File does not exist: " + bpPath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("File does not exist: " + pointPath).toStdString();
    }
}

void CemrgCarpUtils::NormaliseFibreFiles(QString fibresPath, QString outputPath){
    MITK_INFO << "Normalise fibres file";
    std::ifstream ffibres(fibresPath.toStdString());
    std::ofstream fo(outputPath.toStdString());

    int numVect;
    ffibres >> numVect;
    MITK_INFO << ("Number of vectors per line in file: " + QString::number(numVect)).toStdString();

    double x,y,z;
    double norm;
    // prime read
    ffibres >> x;
    ffibres >> y;
    ffibres >> z;
    fo << numVect << std::endl;
    while (!ffibres.eof()) {
        for (int i=0; i<numVect; i++) {
            norm=sqrt(x*x + y*y +z*z);
            if(norm>0){
                fo <<std::fixed << std::setprecision(8) << x/norm << " " << y/norm << " " << z/norm << " ";
            } else{
                fo <<std::fixed << std::setprecision(8) << x << " " << y << " " << z << " ";
            }
            ffibres >> x;
            ffibres >> y;
            ffibres >> z;
        }
        fo << std::endl;
    }
    ffibres.close();
    fo.close();
}

void CemrgCarpUtils::CarpToVtk(QString elemPath, QString ptsPath, QString outputPath, bool saveRegionlabels){
    std::ofstream VTKFile;
    std::ifstream ptsFileRead, elemFileRead;
    short int precision=12;
    short int numColsLookupTable=1;

    VTKFile.open(outputPath.toStdString());
    MITK_INFO << "Writing vtk file header.";
    VTKFile << "# vtk DataFile Version 4.0"<< std::endl;
    VTKFile << "vtk output"<< std::endl;
    VTKFile << "ASCII" << std::endl;
    VTKFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    int nElem, nPts;
    ptsFileRead.open(ptsPath.toStdString());
    ptsFileRead >> nPts;

    MITK_INFO << "Setting geometry - Points";
    VTKFile<<"POINTS "<< nPts <<" float"<<std::endl;
    double x, y, z;
    for (int ix = 0; ix < nPts; ix++) {
        ptsFileRead >> x;
        ptsFileRead >> y;
        ptsFileRead >> z;

        VTKFile<<std::setprecision(precision)<<x<<" "<<y<<" "<<z<<std::endl;
    }
    ptsFileRead.close();

    MITK_INFO << "Setting geometry - Tetrahedral elements";
    elemFileRead.open(elemPath.toStdString());
    elemFileRead >> nElem;
    std::string type;
    int p0, p1, p2, p3, region;
    std::vector<double> regionVector(nElem);
    VTKFile << "CELLS " << nElem << " " << (4+1)*nElem << std::endl;
    for (int ix = 0; ix < nElem; ix++) {
        elemFileRead >> type;
        elemFileRead >> p0;
        elemFileRead >> p1;
        elemFileRead >> p2;
        elemFileRead >> p3;
        elemFileRead >> regionVector[ix];

        VTKFile << "4 " << p0 << " " << p1 << " " << p2 << " " << p3 << std::endl;
    }

    VTKFile <<"CELL_TYPES "<< nElem << " ";
    for (int ix = 0; ix < nElem; ix++) {
        VTKFile << "10"; // type for tetrahedral mesh
        if(((1+ix)%numColsLookupTable) && (ix<(nElem-1))) {
            VTKFile << " ";
        } else{
            VTKFile << std::endl;
        }
    }
    elemFileRead.close();
    VTKFile.close();
    if(saveRegionlabels){
        AppendScalarFieldToVtk(outputPath, "region_labels", "CELL",  regionVector);
    }
}

void CemrgCarpUtils::RectifyFileValues(QString pathToFile, double minVal, double maxVal){
    QFileInfo fi(pathToFile);
    QString copyName = fi.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + fi.baseName() + "_copy." + fi.completeSuffix();
    MITK_INFO << "Copying path name";
    QFile::rename(pathToFile, copyName);

    std::ifstream readInFile(copyName.toStdString());
    std::ofstream writeOutFile(pathToFile.toStdString());

    double valueIn, valueOut;
    int precision = 16;
    int count = 0;
    writeOutFile << std::scientific;
    MITK_INFO << "Rectifying file.";
    while (!readInFile.eof()) {
        valueIn = -1;
        readInFile >> valueIn;
        if(valueIn > maxVal){
            valueOut = maxVal;
        } else if(valueIn < minVal){
            valueOut = minVal;
        } else {
            valueOut = valueIn;
        }
        if(valueIn != -1){
            writeOutFile << std::setprecision(precision) << valueOut << std::endl;
            count++;
        }
    }
    MITK_INFO << ("Finished rectifying file with :" + QString::number(count) + " points.").toStdString();
    readInFile.close();
    writeOutFile.close();

    QFile::remove(copyName);
}

int CemrgCarpUtils::GetTotalFromCarpFile(QString pathToFile, bool totalAtTop){
    int total=-1;
    std::ifstream fi(pathToFile.toStdString());
    if(totalAtTop){
        fi >> total;
    } else{
        int count;
        double value;
        count=0;
        while (!fi.eof()){
            value = -1;
            fi >> value;
            count++;
        }
        total = count;
        total -= (value==-1) ? 1 : 0; // checks if last line in file is empty
    }
    fi.close();
    return total;
}

std::vector<double> CemrgCarpUtils::ReadScalarField(QString pathToFile){
    std::ifstream fi(pathToFile.toStdString());
    int n = CemrgCarpUtils::GetTotalFromCarpFile(pathToFile, false);
    std::vector<double> field(n, 0.0);

    for (int ix = 0; ix < n; ix++) {
        if(fi.eof()){
            MITK_INFO << "File finished prematurely.";
            break;
        }
        fi >> field[ix];
    }
    fi.close();

    return field;
}

void CemrgCarpUtils::AppendScalarFieldToVtk(QString vtkPath, QString fieldName, QString typeData, std::vector<double> field, bool setHeader){
    std::ofstream VTKFile;
    short int precision=12;
    short int numColsLookupTable=1;

    VTKFile.open(vtkPath.toStdString(), std::ios_base::app);
    int fieldSize=field.size();

    if(setHeader){
        MITK_INFO << "Setting POINT_DATA header.";
        VTKFile << typeData.toStdString() << "_DATA " << fieldSize << std::endl;
    }

    MITK_INFO << ("Appending scalar field <<" + fieldName + ">> to VTK file.").toStdString();
    VTKFile << "SCALARS " << fieldName.toStdString() << " FLOAT " << numColsLookupTable<< " " <<std::endl;
    VTKFile << "LOOKUP_TABLE default " << std::endl;

    for (int ix = 0; ix < fieldSize; ix++) {
        VTKFile << std::setprecision(precision) << field.at(ix);
        if(((1+ix)%numColsLookupTable) && (ix<fieldSize-1)){
          VTKFile<<" ";
        } else{
          VTKFile<<std::endl;
        }
    }
    VTKFile.close();
}

void CemrgCarpUtils::AppendVectorFieldToVtk(QString vtkPath, QString fieldName, QString dataType,  std::vector<double> field, bool setHeader){
    std::ofstream VTKFile;
    short int precision=12;
    // short int numColsLookupTable=1;

    VTKFile.open(vtkPath.toStdString(), std::ios_base::app);
    int nElem=field.size()/3;

    if(setHeader){
        MITK_INFO << "Setting CELL_DATA header.";
        VTKFile << dataType.toStdString() << "_DATA " << nElem << std::endl;
    }

    MITK_INFO << ("Appending vector field <<" + fieldName + ">> to VTK file.").toStdString();
    VTKFile << "VECTORS " << fieldName.toStdString() << " FLOAT " << std::endl;

    double x,y,z;
    for (int ix = 0; ix < nElem; ix++) {
        x = field[ix + 0*nElem];
        y = field[ix + 1*nElem];
        z = field[ix + 2*nElem];

        VTKFile << std::setprecision(precision)<<x<<" "<<y<<" "<<z<<std::endl;
    }

    VTKFile.close();
}
