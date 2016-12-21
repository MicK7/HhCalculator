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
#include "vtkSMPTools.h"
#include "vtkSMPThreadLocalObject.h"
#include "vtkDataArray.h"
#include "vtkArrayDispatch.h"
#include "vtkDataArrayAccessor.h"
#include "vtkGenericDataArray.h"
#include "vtkAssume.h"

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
    if (!this->KdTree) {
        vtkPointSet *pointSet = vtkPointSet::SafeDownCast(this->GetDataSet());
        if (!pointSet) {
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
    this->DownFileName = NULL;
    this->TopFileName = NULL;

    //
    this->lineA = vtkSmartPointer<vtkPolyData>::New();
    this->lineB = vtkSmartPointer<vtkPolyData>::New();
    this->NormalsLineA = vtkSmartPointer<vtkDoubleArray>::New();
    this->NormalsLineB = vtkSmartPointer<vtkDoubleArray>::New();

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

void vtkHhCalculator::loadPoints(const char *FileName, vtkSmartPointer<vtkPoints> sortedPoints)
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
    if (!this->DownFileName)
    {
        return 0;
    }
    if (!this->TopFileName)
    {
        return 0;
    }
    if (!vtksys::SystemTools::FileExists(this->DownFileName))
    {
        return 0;
    }
    if (!vtksys::SystemTools::FileExists(this->TopFileName))
    {
        return 0;
    }

    if (this->DownFileName != this->LastLineAFile)
    {
        vtkSmartPointer<vtkPoints> sortedPointsA = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> cellsA = vtkSmartPointer<vtkCellArray>::New();
        vtkIdTypeArray *cellLocationsA = vtkIdTypeArray::New();  // Init a vtkIdTypeArray with 1 component
        //
        this->loadPoints(this->DownFileName, sortedPointsA);
        //
        vtkIdType num_cell = sortedPointsA->GetNumberOfPoints() - 1;
        cellLocationsA->SetNumberOfValues(num_cell * 3);  // NumValue/numComp (==1) - equivalent to SetNumberOfTuples ?

        for (vtkIdType ii = 0; ii < num_cell; ii++)
        {
            // SetComponent instead of SetTuple1
            // because we don't want to check NumberOfComponents each iteration
            cellLocationsA->SetComponent(ii * 3, 0,    2);
            cellLocationsA->SetComponent(ii * 3 + 1, 0,   ii);
            cellLocationsA->SetComponent(ii * 3 + 2, 0, ii + 1);
        }

        cellsA->SetCells(num_cell, cellLocationsA);
        cellLocationsA->Delete();
        //
        this->lineA->SetPoints(sortedPointsA);
        this->lineA->SetLines(cellsA);
        //
        this->NormalsLineA->SetNumberOfComponents(3);
        this->NormalsLineA->SetNumberOfTuples(num_cell);
        //

        for (vtkIdType ii = 0; ii < num_cell; ++ii)
        {
            double pt_cur[3];
            double pt_plus[3];
            double normal[3];
            double norm;
            normal[2] = 0.0;
            //
            sortedPointsA->GetPoint(ii, pt_cur);
            sortedPointsA->GetPoint((ii + 1), pt_plus);
            //
            normal[1] = pt_plus[0] - pt_cur[0];
            normal[0] = - (pt_plus[1] - pt_cur[1]);
            norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
            if (norm < 1.0e-12) {
                norm = 1.0;
            }
            normal[0] /= norm;
            normal[1] /= norm;
            this->NormalsLineA->SetTypedTuple(ii, normal);
        }

        this->LastLineAFile = this->DownFileName;
    }

    if (this->TopFileName != this->LastLineBFile)
    {
        //
        vtkSmartPointer<vtkPoints> sortedPointsB = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> cellsB = vtkSmartPointer<vtkCellArray>::New();
        vtkIdTypeArray *cellLocationsB = vtkIdTypeArray::New();

        this->loadPoints(this->TopFileName, sortedPointsB);
        vtkIdType num_cell = sortedPointsB->GetNumberOfPoints() - 1;
        //
        cellLocationsB->SetNumberOfValues(num_cell * 3);
        for (vtkIdType ii = 0; ii < num_cell; ii++)
        {
            // SetComponent instead of SetTuple1
            // because we don't want to check NumberOfComponents each iteration
            cellLocationsB->SetComponent(ii * 3, 0,    2);
            cellLocationsB->SetComponent(ii * 3 + 1, 0,   ii);
            cellLocationsB->SetComponent(ii * 3 + 2, 0, ii + 1);
        }
        cellsB->SetCells(num_cell, cellLocationsB);
        cellLocationsB->Delete();
        //
        this->lineB->SetPoints(sortedPointsB);
        this->lineB->SetLines(cellsB);
        //
        this->NormalsLineB->SetNumberOfComponents(3);
        this->NormalsLineB->SetNumberOfTuples(num_cell);
        for (vtkIdType ii = 0; ii < num_cell; ++ii)
        {
            double pt_cur[3];
            double pt_plus[3];
            double normal[3];
            double norm;
            normal[2] = 0.0;
            //
            sortedPointsB->GetPoint(ii, pt_cur);
            sortedPointsB->GetPoint(ii + 1, pt_plus);
            //
            normal[1] = - (pt_plus[0] - pt_cur[0]);
            normal[0] = + (pt_plus[1] - pt_cur[1]);
            norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
            if (norm < 1.0e-12)
            {
                norm = 1.0;
            }
            normal[0] /= norm;
            normal[1] /= norm;
            this->NormalsLineB->SetTypedTuple(ii, normal);
        }
        this->LastLineBFile = this->TopFileName;
    }

    return 1;
}

//----------------------------------------------------------------------------
inline double computeCurveIndex(const double *pt_p,
                                const double *pt_m,
                                const double *pt_r)
{
    double t = 0.0;
    double diff[4];
    diff[0] = pt_p[0] - pt_m[0];
    diff[1] = pt_p[1] - pt_m[1];
    diff[2] = pt_r[0] - pt_m[0];
    diff[3] = pt_r[1] - pt_m[1];
    t = (diff[0] * diff[2] + diff[1] * diff[3]) / (diff[0] * diff[0] + diff[1] * diff[1]);
    return t;
}

//----------------------------------------------------------------------------
template <typename PointArray, typename DistanceArray>
class computeHhFunctor
{
public:
    computeHhFunctor(double Scaling, double Tolerance,
            PointArray *points,
            DistanceArray *distanceA,
            DistanceArray *distanceB,
            DistanceArray *dOverD,
            vtkPolyData *lineA,
            vtkPolyData *lineB,
            vtkDoubleArray *NormalsLineA,
            vtkDoubleArray *NormalsLineB) :
        Scaling(Scaling), Tolerance(Tolerance), points(points), distanceA(distanceA), distanceB(distanceB), dOverD(dOverD),
        lineA(lineA), lineB(lineB), NormalsLineA(NormalsLineA), NormalsLineB(NormalsLineB)
    {
    }

    void Initialize()
    {
        ThreadLocalWorkSpace &ws = this->tlws.Local();
        ws.pointLocatorA = vtkSmartPointer<vtkKdTree2dPointLocator>::New();
        ws.pointLocatorB = vtkSmartPointer<vtkKdTree2dPointLocator>::New();
        ws.lineA = vtkSmartPointer<vtkPolyData>::New();
        ws.lineA->DeepCopy(this->lineA);
        ws.lineB = vtkSmartPointer<vtkPolyData>::New();
        ws.lineB->DeepCopy(this->lineB);
        //
        ws.pointLocatorA->SetDataSet(ws.lineA);
        ws.pointLocatorA->Build2dLocator();
        //
        ws.pointLocatorB->SetDataSet(ws.lineB);
        ws.pointLocatorB->Build2dLocator();
    }

    void operator()(vtkIdType begin, vtkIdType end)
    {
        //
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
        //
        // This allows the compiler to optimize for the AOS array stride.
        VTK_ASSUME(this->points->GetNumberOfComponents() == 3);
        VTK_ASSUME(this->distanceA->GetNumberOfComponents() == 1);
        VTK_ASSUME(this->distanceB->GetNumberOfComponents() == 1);
        VTK_ASSUME(this->dOverD->GetNumberOfComponents() == 1);

        vtkDataArrayAccessor<PointArray> pt(points);
        vtkDataArrayAccessor<DistanceArray> dA(this->distanceA);
        vtkDataArrayAccessor<DistanceArray> dB(this->distanceB);
        vtkDataArrayAccessor<DistanceArray> dD(this->dOverD);

        // Get pointLocatorA and pointLocatorB
        ThreadLocalWorkSpace &ws = tlws.Local();
        vtkKdTree2dPointLocator *pointLocatorA = ws.pointLocatorA.GetPointer();
        vtkKdTree2dPointLocator *pointLocatorB = ws.pointLocatorB.GetPointer();
        //
        vtkIdType Seg0PointPId;
        vtkIdType Seg0PointMId;
        vtkIdType Seg1PointPId;
        vtkIdType Seg1PointMId;

        double tA_seg0, tA_seg1;
        double tB_seg0, tB_seg1;
        int nPtsA;
        int nPtsB;
        const double Ldim = this->Scaling;
        const double Tolerance = this->Tolerance;
        //
        nPtsA = this->lineA->GetNumberOfPoints();
        nPtsB = this->lineB->GetNumberOfPoints();

        for (vtkIdType tupleIdx = begin; tupleIdx < end; ++tupleIdx)
        {
            double currentPoint[3];
            double currentXRPoint[3];
            currentXRPoint[0] = currentXRPoint[1] = currentXRPoint[2] = 0.0;
            // Get Current XR point
            currentPoint[0] = pt.Get(tupleIdx, 0);
            currentPoint[1] = pt.Get(tupleIdx, 1);
            currentPoint[2] = pt.Get(tupleIdx, 2);
            currentXRPoint[0] = Ldim * currentPoint[0];
            currentXRPoint[1] = Ldim * std::sqrt(currentPoint[1]*currentPoint[1]+currentPoint[2]*currentPoint[2]);
            //
            vtkIdType closestPointAId = pointLocatorA->FindClosestPoint(currentXRPoint);
            vtkIdType closestPointBId = pointLocatorB->FindClosestPoint(currentXRPoint);
            //---------
            // lineA
            if (closestPointAId == 0)
            {
                Seg0PointMId = closestPointAId;
                Seg1PointMId = closestPointAId;
                Seg0PointPId = closestPointAId + 1;
                Seg1PointPId = closestPointAId + 1;
            }
            else if (closestPointAId == (nPtsA - 1))
            {
                Seg0PointMId = closestPointAId - 1;
                Seg1PointMId = closestPointAId - 1;
                Seg0PointPId = closestPointAId;
                Seg1PointPId = closestPointAId;
            }
            else
            {
                Seg0PointMId = closestPointAId - 1;
                Seg1PointMId = closestPointAId;
                Seg0PointPId = closestPointAId;
                Seg1PointPId = closestPointAId + 1;
            }

            //
            this->lineA->GetPoint(Seg0PointMId, pt_m_seg0);
            this->lineA->GetPoint(Seg0PointPId, pt_p_seg0);
            this->lineA->GetPoint(Seg1PointMId, pt_m_seg1);
            this->lineA->GetPoint(Seg1PointPId, pt_p_seg1);
            this->lineA->GetPoint(closestPointAId, closestPointA);
            //
            tA_seg0 = computeCurveIndex(pt_p_seg0, pt_m_seg0, currentXRPoint);
            tA_seg1 = computeCurveIndex(pt_p_seg1, pt_m_seg1, currentXRPoint);
            //
            this->NormalsLineA->GetTypedTuple(Seg0PointMId, norms_seg0);
            this->NormalsLineA->GetTypedTuple(Seg1PointMId, norms_seg1);

            // For segment 0
            distance_proj0 = norms_seg0[0] * (currentXRPoint[0] - closestPointA[0])
                             + norms_seg0[1] * (currentXRPoint[1] - closestPointA[1]);

            if (tA_seg0 >= 1.0 || tA_seg0 < 0.0)
            {
                distance_seg0 = std::sqrt((currentXRPoint[0] - closestPointA[0]) * (currentXRPoint[0] - closestPointA[0])
                                          + (currentXRPoint[1] - closestPointA[1]) * (currentXRPoint[1] - closestPointA[1]));
            }
            else
            {
                distance_seg0 = std::fabs(distance_proj0);
            }
            // For segment 1
            distance_proj1 = norms_seg1[0] * (currentXRPoint[0] - closestPointA[0])
                             + norms_seg1[1] * (currentXRPoint[1] - closestPointA[1]);


            if (tA_seg1 > 1.0 || tA_seg1 <= 0.0)
            {
                distance_seg1 = std::sqrt((currentXRPoint[0] - closestPointA[0]) * (currentXRPoint[0] - closestPointA[0])
                                          + (currentXRPoint[1] - closestPointA[1]) * (currentXRPoint[1] - closestPointA[1]));
            }
            else
            {
                distance_seg1 = std::fabs(distance_proj1);
            }
            //
            if (distance_seg0 < distance_seg1)
            {
                signed_distanceA = distance_seg0;
                if (distance_seg0 < Tolerance)
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
                if (distance_seg1 < Tolerance)
                {
                    signed_distanceA = 0.0;
                }
                if (distance_proj1 < 0.0)
                {
                    signed_distanceA *= -1.0;
                }
            }
            //
            dA.Set(tupleIdx, 0, signed_distanceA);
            //
            //---------
            // lineB
            if (closestPointBId == 0)
            {
                Seg0PointMId = closestPointBId;
                Seg1PointMId = closestPointBId;
                Seg0PointPId = closestPointBId + 1;
                Seg1PointPId = closestPointBId + 1;
            }
            else if (closestPointBId == (nPtsB - 1))
            {
                Seg0PointMId = closestPointBId - 1;
                Seg1PointMId = closestPointBId - 1;
                Seg0PointPId = closestPointBId;
                Seg1PointPId = closestPointBId;
            }
            else
            {
                Seg0PointMId = closestPointBId - 1;
                Seg1PointMId = closestPointBId;
                Seg0PointPId = closestPointBId;
                Seg1PointPId = closestPointBId + 1;
            }
            //
            this->lineB->GetPoint(Seg0PointMId, pt_m_seg0);
            this->lineB->GetPoint(Seg1PointMId, pt_m_seg1);
            this->lineB->GetPoint(Seg0PointPId, pt_p_seg0);
            this->lineB->GetPoint(Seg1PointPId, pt_p_seg1);
            this->lineB->GetPoint(closestPointBId, closestPointB);
            //
            tB_seg0 = computeCurveIndex(pt_p_seg0, pt_m_seg0, currentXRPoint);
            tB_seg1 = computeCurveIndex(pt_p_seg1, pt_m_seg1, currentXRPoint);
            //
            this->NormalsLineB->GetTypedTuple(Seg0PointMId, norms_seg0);
            this->NormalsLineB->GetTypedTuple(Seg1PointMId, norms_seg1);
            // For segment 0
            distance_proj0 = norms_seg0[0] * (currentXRPoint[0] - closestPointB[0])
                             + norms_seg0[1] * (currentXRPoint[1] - closestPointB[1]);
            //
            if (tB_seg0 >= 1.0 || tB_seg0 < 0.0)
            {
                distance_seg0 = std::sqrt((currentXRPoint[0] - closestPointB[0]) * (currentXRPoint[0] - closestPointB[0])
                                          + (currentXRPoint[1] - closestPointB[1]) * (currentXRPoint[1] - closestPointB[1]));
            }
            else
            {
                distance_seg0 = std::fabs(distance_proj0);
            }
            //
            // For segment 1
            distance_proj1 = norms_seg1[0] * (currentXRPoint[0] - closestPointB[0])
                             + norms_seg1[1] * (currentXRPoint[1] - closestPointB[1]);

            if (tB_seg1 > 1.0 || tB_seg1 <= 0.0)
            {
                distance_seg1 = std::sqrt((currentXRPoint[0] - closestPointB[0]) * (currentXRPoint[0] - closestPointB[0])
                                          + (currentXRPoint[1] - closestPointB[1]) * (currentXRPoint[1] - closestPointB[1]));
            }
            else
            {
                distance_seg1 = std::fabs(distance_proj1);
            }
            if (distance_seg0 < distance_seg1)
            {
                signed_distanceB = distance_seg0;
                if (distance_seg0 < Tolerance)
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
                if (distance_seg1 < Tolerance)
                {
                    signed_distanceB = 0.0;
                }
                if (distance_proj1 < 0.0)
                {
                    signed_distanceB *= -1.0;
                }
            }
            //
            dB.Set(tupleIdx, 0, signed_distanceB);
            //
            // compute dOverD
            double dist = (signed_distanceA + signed_distanceB);
            // .. Buggy
            if (dist == 0)
            {
                dist = 1.0e-3;
            }
            dist = signed_distanceA / dist;
            dD.Set(tupleIdx, 0, dist);
        }
    }

    void Reduce()
    {
    }

private:
    struct ThreadLocalWorkSpace
    {
        vtkSmartPointer<vtkKdTree2dPointLocator> pointLocatorA;
        vtkSmartPointer<vtkKdTree2dPointLocator> pointLocatorB;
        vtkSmartPointer<vtkPolyData> lineA;
        vtkSmartPointer<vtkPolyData> lineB;
    };

    typedef vtkSMPThreadLocal<ThreadLocalWorkSpace> TLS_t;
    TLS_t tlws;

    double Scaling, Tolerance;
    PointArray *points;
    DistanceArray *distanceA;
    DistanceArray *distanceB;
    DistanceArray *dOverD;
    vtkDoubleArray *NormalsLineA;
    vtkDoubleArray *NormalsLineB;
    vtkPolyData *lineA;
    vtkPolyData *lineB;

};

//----------------------------------------------------------------------------
template <typename DistanceArray>
class HhWorker
{
private:
    DistanceArray *distanceA;
    DistanceArray *distanceB;
    DistanceArray *dOverD;
    double Scaling, Tolerance;
    vtkPolyData *lineA;
    vtkPolyData *lineB;
    vtkDoubleArray *NormalsLineA;
    vtkDoubleArray *NormalsLineB;

public:
    HhWorker (DistanceArray *distanceA,
              DistanceArray *distanceB,
              DistanceArray *dOverD,
              const double Scaling, const double Tolerance,
              vtkPolyData *lineA,
              vtkPolyData *lineB,
              vtkDoubleArray *NormalsLineA,
              vtkDoubleArray *NormalsLineB) :
        distanceA(distanceA), distanceB(distanceB), dOverD(dOverD), Scaling(Scaling), Tolerance(Tolerance),
        lineA(lineA), lineB(lineB), NormalsLineA(NormalsLineA), NormalsLineB(NormalsLineB)
    {
    }

    // The worker accepts VTK array objects, not raw memory buffers.
    template <typename PointArray>
    void operator()(PointArray *points)
    {
        // This allows the compiler to optimize for the AOS array stride.
        VTK_ASSUME(points->GetNumberOfComponents() == 3);

        vtkIdType numPoints = points->GetNumberOfTuples();
        // These allow this single worker function to be used with both
        // the vtkDataArray 'double' API and the more efficient
        // vtkGenericDataArray APIs, depending on the template parameters:
        computeHhFunctor<PointArray, DistanceArray>  functorHh(this->Scaling, this->Tolerance,
                points,
                this->distanceA,
                this->distanceB,
                this->dOverD,
                this->lineA,
                this->lineB,
                this->NormalsLineA,
                this->NormalsLineB);
        vtkSMPTools::For(0, numPoints, functorHh);
    }
};

//----------------------------------------------------------------------------
int vtkHhCalculator::RequestData(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkDataSet *in = vtkDataSet::SafeDownCast(
                         inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataSet *out = vtkDataSet::SafeDownCast(
                          outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Copy input to get a start point
    out->ShallowCopy(in);

    if (this->UpdateLocators() != VTK_OK)
    {
        return 0;
    }

    vtkPointSet *inPtSet = static_cast<vtkPointSet *>(in);
    if (! inPtSet)
    {
        std::cerr << "Non vtkPointSet are not supported right now" << std::endl;
        return 0;
    }

    // Allocate storage for d_over_D
    vtkIdType const nPts = in->GetNumberOfPoints();
    vtkDataArray *pts = inPtSet->GetPoints()->GetData();

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

    HhWorker<vtkDoubleArray> worker(distanceA.GetPointer(),
                                    distanceB.GetPointer(),
                                    dOverD.GetPointer(),
                                    this->Scaling, this->Tolerance,
                                    lineA.GetPointer(), lineB.GetPointer(),
                                    NormalsLineA.GetPointer(), NormalsLineB.GetPointer());

    // Define our dispatcher. We’ll only consider float/double arrays.
    // These combinations will use a 'fast-path'
    // implementation generated by the dispatcher:
    typedef vtkArrayDispatch::DispatchByValueType
    <
    vtkArrayDispatch::Reals // ValueTypes allowed by first array points
    > Dispatcher;

    // Execute the dispatcher:
    if (!Dispatcher::Execute(pts, worker))
    {
        // If Execute() fails, it means the dispatch failed due to an
        // unsupported array type. In this case, it’s likely that the points
        // array is using an integral type. This is an uncommon case, so we won’t
        // generate a fast path for these, but instead call an instantiation of
        // Through the use of vtkDataArrayAccessor, this falls back to using the
        // vtkDataArray double API:
        worker(pts);
    }

    //----
    out->GetPointData()->AddArray(dOverD);
    out->GetPointData()->SetActiveAttribute("dOverD", vtkDataSetAttributes::SCALARS);

    return 1;
}
