/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkHhCalculator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//  Copyright 2012-2016 Mickael Philit.

#include "vtkHhCalculator.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include <fstream>
#include <string>
#include <iostream>

#include <vtksys/SystemTools.hxx>

#include "vtkObjectFactory.h"
#include "vtkDataSetAttributes.h"
#include "vtkPointData.h"
#include "vtkKdTree.h"
//
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkKdTree2dPointLocator);

//----------------------------------------------------------------------------
vtkKdTree2dPointLocator::vtkKdTree2dPointLocator()
{

}

//----------------------------------------------------------------------------
vtkKdTree2dPointLocator::~vtkKdTree2dPointLocator()
{
}

//----------------------------------------------------------------------------
void vtkKdTree2dPointLocator::Build2dLocator()
{
  if(!this->KdTree)
    {
    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(this->GetDataSet());
    if(!pointSet)
      {
      vtkErrorMacro("vtkKdTreePointLocator requires a PointSet to build locator.");
      return;
      }
    this->KdTree = vtkKdTree::New();
    this->KdTree->OmitZPartitioning();    
    this->KdTree->BuildLocatorFromPoints(pointSet);
    this->KdTree->GetBounds(this->Bounds);
    this->Modified();
    }
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkHhCalculator);

//----------------------------------------------------------------------------
vtkHhCalculator::vtkHhCalculator()
{
    this->HubFileName = NULL;
    this->ShroudFileName = NULL;

    //
    this->lineA = vtkSmartPointer<vtkPolyData>::New();
    this->lineB = vtkSmartPointer<vtkPolyData>::New();
    this->NormalsLineA = vtkSmartPointer<vtkDoubleArray>::New();
    this->NormalsLineB = vtkSmartPointer<vtkDoubleArray>::New();
    this->pointLocatorA = vtkSmartPointer<vtkKdTree2dPointLocator>::New();
    this->pointLocatorB = vtkSmartPointer<vtkKdTree2dPointLocator>::New();
    //
    this->LastLineAFile = NULL;
    this->LastLineBFile = NULL;

    this->Tolerance = 0.0005;
    this->Scaling = 1000.0;
    //this->DebugOn();
}

//----------------------------------------------------------------------------
vtkHhCalculator::~vtkHhCalculator()
{

}

//----------------------------------------------------------------------------
void vtkHhCalculator::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "vtkHhCalculator" << std::endl;
}

void vtkHhCalculator::loadPoints(const char* FileName, vtkSmartPointer<vtkPoints> sortedPoints)
{
    std::string Name;
    std::string Type;
    int NumberOfPoints;
    std::fstream fs;
    fs.open(FileName, std::fstream::in);
    fs >> Name;
    fs >> Type;
    if (Type != "ZR")
    {
        sortedPoints->SetNumberOfPoints(0);
        fs.close();
        return;
    }
    fs >> NumberOfPoints;
    sortedPoints->SetDataTypeToDouble();
    sortedPoints->SetNumberOfPoints(NumberOfPoints);
    for (int ii=0; ii <NumberOfPoints; ii++)
    {
        double a, b;
        fs >> a >> b;
        sortedPoints->SetPoint(ii, a, b, 0.0);
    }
    fs.close();
}

//----------------------------------------------------------------------------
int vtkHhCalculator::UpdateLocators()
{

    // check if filename are not NULL and FileExists
    if (!this->HubFileName)
    {
        return 0;
    }
    if (!this->ShroudFileName)
    {
        return 0;
    }
    if (!vtksys::SystemTools::FileExists(this->HubFileName))
    {
        return 0;
    }
    if (!vtksys::SystemTools::FileExists(this->ShroudFileName))
    {
        return 0;
    }

    if (this->HubFileName != this->LastLineAFile)
    {
        vtkSmartPointer<vtkPoints> sortedPointsA = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> cellsA = vtkSmartPointer<vtkCellArray>::New();
        vtkIdTypeArray * cellLocationsA = vtkIdTypeArray::New();
        //
        this->loadPoints(this->HubFileName, sortedPointsA);
        //
        cellLocationsA->SetNumberOfValues((sortedPointsA->GetNumberOfPoints()-1)*3);
        for (vtkIdType ii=0; ii<sortedPointsA->GetNumberOfPoints()-1; ii++)
        {
            cellLocationsA->SetTuple1(ii*3,2);
            cellLocationsA->SetTuple1(ii*3+1,ii);
            cellLocationsA->SetTuple1(ii*3+2,ii+1);
        }
        cellsA->SetCells(sortedPointsA->GetNumberOfPoints()-1,cellLocationsA);
        cellLocationsA->Delete();
        //
        this->lineA->SetPoints(sortedPointsA);
        this->lineA->SetLines(cellsA);
        //
        this->NormalsLineA->SetNumberOfComponents(3);
        this->NormalsLineA->SetNumberOfTuples(sortedPointsA->GetNumberOfPoints()-1);
        //

        for (vtkIdType ii=0; ii<sortedPointsA->GetNumberOfPoints()-1; ++ii)
        {
            double pt_cur[3];
            double pt_plus[3];
            double normal[3];
            double norm;
            normal[2] = 0.0;
            //
            sortedPointsA->GetPoint(ii, pt_cur);
            sortedPointsA->GetPoint((ii+1),pt_plus);
            //
            normal[1] = pt_plus[0]-pt_cur[0];
            normal[0] = -(pt_plus[1]-pt_cur[1]);
            norm = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
            if (norm < 1.0e-12)
            {
                norm = 1.0;
            }
            normal[0] /= norm;
            normal[1] /= norm;
            this->NormalsLineA->SetTuple(ii,normal);
        }

        this->pointLocatorA->SetDataSet(this->lineA);
        this->pointLocatorA->Build2dLocator();
        this->LastLineAFile = this->HubFileName;
    }

    if (this->ShroudFileName != this->LastLineBFile)
    {
        //
        vtkSmartPointer<vtkPoints> sortedPointsB = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> cellsB = vtkSmartPointer<vtkCellArray>::New();
        vtkIdTypeArray* cellLocationsB = vtkIdTypeArray::New();

        this->loadPoints(this->ShroudFileName, sortedPointsB);
        //
        cellLocationsB->SetNumberOfValues((sortedPointsB->GetNumberOfPoints()-1)*3);
        for (vtkIdType ii=0; ii<sortedPointsB->GetNumberOfPoints()-1; ii++)
        {
            cellLocationsB->SetTuple1(ii*3,2);
            cellLocationsB->SetTuple1(ii*3+1,ii);
            cellLocationsB->SetTuple1(ii*3+2,ii+1);
        }
        cellsB->SetCells(sortedPointsB->GetNumberOfPoints()-1,cellLocationsB);
        cellLocationsB->Delete();
        //
        this->lineB->SetPoints(sortedPointsB);
        this->lineB->SetLines(cellsB);
        //
        this->NormalsLineB->SetNumberOfComponents(3);
        this->NormalsLineB->SetNumberOfTuples(sortedPointsB->GetNumberOfPoints()-1);
        for (vtkIdType ii=0; ii<sortedPointsB->GetNumberOfPoints()-1; ++ii)
        {
            double pt_cur[3];
            double pt_plus[3];
            double normal[3];
            double norm;
            normal[2] = 0.0;
            //
            sortedPointsB->GetPoint(ii, pt_cur);
            sortedPointsB->GetPoint(ii+1, pt_plus);
            //
            normal[1] = -(pt_plus[0]-pt_cur[0]);
            normal[0] = +(pt_plus[1]-pt_cur[1]);
            norm = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
            if (norm < 1.0e-12)
            {
                norm = 1.0;
            }
            normal[0] /= norm;
            normal[1] /= norm;
            this->NormalsLineB->SetTuple(ii,normal);
        }

        //this->pointLocatorB->KdTree->OmitZPartitioning();
        this->pointLocatorB->SetDataSet(this->lineB);
        this->pointLocatorB->Build2dLocator();
        this->LastLineBFile = this->ShroudFileName;
    }

    return 1;
}

//----------------------------------------------------------------------------
inline double computeCurveIndex(const double *pt_p,
                                const double* pt_m,
                                const double *pt_r)
{
    double t = 0.0;
    double diff[4];
    diff[0] = pt_p[0]-pt_m[0];
    diff[1] = pt_p[1]-pt_m[1];
    diff[2] = pt_r[0]-pt_m[0];
    diff[3] = pt_r[1]-pt_m[1];
    t = (diff[0]*diff[2] + diff[1]*diff[3])/(diff[0]*diff[0]+diff[1]*diff[1]);
    return t;
}

//----------------------------------------------------------------------------
int vtkHhCalculator::RequestData(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkDataSet* in = vtkDataSet::SafeDownCast(
                inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataSet* out = vtkDataSet::SafeDownCast(
                outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Copy input to get a start point
    out->ShallowCopy(in);

    if ( this->UpdateLocators() != VTK_OK )
    {
        return 0;
    }

    // Allocate storage for d_over_D
    vtkIdType const nPts = in->GetNumberOfPoints();
    vtkSmartPointer<vtkDoubleArray> dOverD =
            vtkSmartPointer<vtkDoubleArray>::New();
    dOverD->SetName("dOverD");
    dOverD->SetNumberOfValues(nPts);
    // Allocate storage for distance to line A and B
    vtkSmartPointer<vtkDoubleArray> distanceA =
            vtkSmartPointer<vtkDoubleArray>::New();
    distanceA->SetName("distanceA");
    distanceA->SetNumberOfValues(nPts);

    vtkSmartPointer<vtkDoubleArray> distanceB =
            vtkSmartPointer<vtkDoubleArray>::New();
    distanceB->SetName("distanceB");
    distanceB->SetNumberOfValues(nPts);

    //
    const double Ldim = this->Scaling;
    vtkIdType Seg0PointPId;
    vtkIdType Seg0PointMId;
    vtkIdType Seg1PointPId;
    vtkIdType Seg1PointMId;

    double tA_seg0, tA_seg1;
    double tB_seg0, tB_seg1;
    int nPtsA;
    int nPtsB;
    //
    nPtsA = this->lineA->GetNumberOfPoints();
    nPtsB = this->lineB->GetNumberOfPoints();
    //

    // Compute distanceA and distanceB --> A == hub, B== shroud
    for(int indexPts = 0; indexPts < nPts; indexPts++)
    {
        double signed_distanceA;
        double signed_distanceB;
        double distance_seg0;
        double distance_seg1;
        double distance_proj0;
        double distance_proj1;
        double closestPointA[3];
        double closestPointB[3];
        double pt_p_seg0[3];
        double pt_m_seg0[3];
        double pt_p_seg1[3];
        double pt_m_seg1[3];
        double norms_seg0[3];
        double norms_seg1[3];
        double currentPoint[3];
        double currentXRPoint[3];
        currentXRPoint[0] = currentXRPoint[1] = currentXRPoint[2] = 0.0;

        //
        in->GetPoint(indexPts, currentPoint);
        currentXRPoint[0] = currentPoint[0]*Ldim;
        currentXRPoint[1] = sqrt(currentPoint[1]*currentPoint[1] + currentPoint[2]*currentPoint[2])*Ldim;
        //
        vtkIdType closestPointAId = pointLocatorA->FindClosestPoint(currentXRPoint);
        vtkIdType closestPointBId = pointLocatorB->FindClosestPoint(currentXRPoint);
        //---------
        // lineA
        if (closestPointAId == 0)
        {
            Seg0PointMId = closestPointAId;
            Seg1PointMId = closestPointAId;
            Seg0PointPId = closestPointAId+1;
            Seg1PointPId = closestPointAId+1;
        }
        else if (closestPointAId == (nPtsA-1))
        {
            Seg0PointMId = closestPointAId-1;
            Seg1PointMId = closestPointAId-1;
            Seg0PointPId = closestPointAId;
            Seg1PointPId = closestPointAId;
        }
        else
        {
            Seg0PointMId = closestPointAId-1;
            Seg1PointMId = closestPointAId;
            Seg0PointPId = closestPointAId;
            Seg1PointPId = closestPointAId+1;
        }
        //
        this->lineA->GetPoint(Seg0PointMId, pt_m_seg0);
        this->lineA->GetPoint(Seg0PointPId, pt_p_seg0);
        this->lineA->GetPoint(Seg1PointMId, pt_m_seg1);
        this->lineA->GetPoint(Seg1PointPId, pt_p_seg1);
        this->lineA->GetPoint(closestPointAId, closestPointA);
        //
        tA_seg0 = computeCurveIndex(pt_p_seg0,pt_m_seg0, currentXRPoint);
        tA_seg1 = computeCurveIndex(pt_p_seg1,pt_m_seg1, currentXRPoint);
        //
        this->NormalsLineA->GetTuple(Seg0PointMId, norms_seg0);
        this->NormalsLineA->GetTuple(Seg1PointMId, norms_seg1);

        // For segment 0
        distance_proj0 = norms_seg0[0]*(currentXRPoint[0]-closestPointA[0])
                       + norms_seg0[1]*(currentXRPoint[1]-closestPointA[1]);

        if (tA_seg0 >= 1.0 || tA_seg0 < 0.0)
        {
            distance_seg0 = sqrt( (currentXRPoint[0]-closestPointA[0])*(currentXRPoint[0]-closestPointA[0])
                                 +(currentXRPoint[1]-closestPointA[1])*(currentXRPoint[1]-closestPointA[1]));
        }
        else
        {
            distance_seg0 = fabs(distance_proj0);
        }
        // For segment 1
        distance_proj1 = norms_seg1[0]*(currentXRPoint[0]-closestPointA[0])
                       + norms_seg1[1]*(currentXRPoint[1]-closestPointA[1]);


        if (tA_seg1 > 1.0 || tA_seg1 <= 0.0)
        {
            distance_seg1 = sqrt((currentXRPoint[0]-closestPointA[0])*(currentXRPoint[0]-closestPointA[0])
                                +(currentXRPoint[1]-closestPointA[1])*(currentXRPoint[1]-closestPointA[1]));
        }
        else
        {
            distance_seg1 = fabs(distance_proj1);
        }
        //
        if (distance_seg0 < distance_seg1)
        {
            signed_distanceA = distance_seg0;
            if (distance_seg0 < this->Tolerance)
            {
                signed_distanceA = 0.0;
            }
            if (distance_proj0 < 0.0)
            {
                signed_distanceA *= -1.0;
            }
        }
        else
        {
            signed_distanceA = distance_seg1;
            if (distance_seg1 < this->Tolerance)
            {
                signed_distanceA = 0.0;
            }
            if (distance_proj1 < 0.0)
            {
                signed_distanceA *= -1.0;
            }
        }

        //
        distanceA->SetTuple1(indexPts, signed_distanceA);
        //
        //---------
        // lineB
        if (closestPointBId == 0)
        {
            Seg0PointMId = closestPointBId;
            Seg1PointMId = closestPointBId;
            Seg0PointPId = closestPointBId+1;
            Seg1PointPId = closestPointBId+1;
        }
        else if (closestPointBId == (nPtsB-1))
        {
            Seg0PointMId = closestPointBId-1;
            Seg1PointMId = closestPointBId-1;
            Seg0PointPId = closestPointBId;
            Seg1PointPId = closestPointBId;
        }
        else
        {
            Seg0PointMId = closestPointBId-1;
            Seg1PointMId = closestPointBId;
            Seg0PointPId = closestPointBId;
            Seg1PointPId = closestPointBId+1;
        }
        //
        this->lineB->GetPoint(Seg0PointMId, pt_m_seg0);
        this->lineB->GetPoint(Seg1PointMId, pt_m_seg1);
        this->lineB->GetPoint(Seg0PointPId, pt_p_seg0);
        this->lineB->GetPoint(Seg1PointPId, pt_p_seg1);
        this->lineB->GetPoint(closestPointBId, closestPointB);
        //
        tB_seg0 = computeCurveIndex(pt_p_seg0,pt_m_seg0,currentXRPoint);
        tB_seg1 = computeCurveIndex(pt_p_seg1,pt_m_seg1,currentXRPoint);
        //
        this->NormalsLineB->GetTuple(Seg0PointMId, norms_seg0);
        this->NormalsLineB->GetTuple(Seg1PointMId, norms_seg1);
        // For segment 0
        distance_proj0 = norms_seg0[0]*(currentXRPoint[0]-closestPointB[0])
                       + norms_seg0[1]*(currentXRPoint[1]-closestPointB[1]);
        //
        if (tB_seg0 >= 1.0 || tB_seg0 < 0.0)
        {
            distance_seg0 = sqrt( (currentXRPoint[0]-closestPointB[0])*(currentXRPoint[0]-closestPointB[0])
                    +(currentXRPoint[1]-closestPointB[1])*(currentXRPoint[1]-closestPointB[1]) );
        }
        else
        {
            distance_seg0 = fabs(distance_proj0);
        }
        //
        // For segment 1
        distance_proj1 = norms_seg1[0]*(currentXRPoint[0]-closestPointB[0])
                       + norms_seg1[1]*(currentXRPoint[1]-closestPointB[1]);

        if (tB_seg1 > 1.0 || tB_seg1 <= 0.0)
        {
            distance_seg1 = sqrt( (currentXRPoint[0]-closestPointB[0])*(currentXRPoint[0]-closestPointB[0])
                    +(currentXRPoint[1]-closestPointB[1])*(currentXRPoint[1]-closestPointB[1]) );
        }
        else
        {
            distance_seg1 = fabs(distance_proj1);
        }
        if (distance_seg0 < distance_seg1)
        {
            signed_distanceB = distance_seg0;
            if (distance_seg0 < this->Tolerance)
            {
                signed_distanceB = 0.0;
            }
            if (distance_proj0 < 0.0)
            {
                signed_distanceB *= -1.0;
            }
        }
        else
        {
            signed_distanceB = distance_seg1;
            if (distance_seg1 < this->Tolerance)
            {
                signed_distanceB = 0.0;
            }
            if (distance_proj1 < 0.0)
            {
                signed_distanceB *= -1.0;
            }
        }
        //
        distanceB->SetTuple1(indexPts, signed_distanceB);
        //
    }
    // Compute dOverD
    for(int i = 0; i < in->GetNumberOfPoints(); i++)
    {
        double dist;
        double dA;
        double dB;
        //
        dA = distanceA->GetTuple1(i);
        dB = distanceB->GetTuple1(i);
        dist = (dA + dB);
        // Buggy
        if (dist  == 0.0)
        {
            dist = 1.0e-3;
        }
        dist = dA/dist;
        dOverD->SetTuple1(i, dist);
    }
    //----
    out->GetPointData()->AddArray(dOverD);
    out->GetPointData()->SetActiveAttribute("dOverD", vtkDataSetAttributes::SCALARS);

    return 1;
}
