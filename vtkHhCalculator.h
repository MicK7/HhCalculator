/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkHhCalculator.h

  Copyright (c) Ken Martin, Will Schrodeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
// Copyright 2012-2016 Mickael Philit

/**
 * @class   vtkHhCalculator
 *
 * vtkHhCalculator parametrize a volume with a distance ratio
 * to axisymmetrical curves
 * 2D curves in meridional plane are defined through a point based file
 *
 * @warning
 *   ...
 *
 * @par Thanks:
 * Thanks to .
 */

#ifndef __vtkHhCalculator_h
#define __vtkHhCalculator_h

#include "vtkSetGet.h"

#include "vtkDataSetAlgorithm.h"
#include "vtkDataSet.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkKdTreePointLocator.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
//
// Helper class
class vtkKdTree2dPointLocator : public vtkKdTreePointLocator
{
public:
    static vtkKdTree2dPointLocator* New();
    vtkTypeMacro(vtkKdTree2dPointLocator, vtkKdTreePointLocator);

    void Build2dLocator();

protected:

    vtkKdTree2dPointLocator();
    ~vtkKdTree2dPointLocator();

private:
    vtkKdTree2dPointLocator(const vtkKdTree2dPointLocator&);
    bool operator=(const vtkKdTree2dPointLocator&);
};

//
class VTK_EXPORT vtkHhCalculator : public vtkDataSetAlgorithm
{
public:
    static vtkHhCalculator* New();
    vtkTypeMacro(vtkHhCalculator, vtkDataSetAlgorithm);
    void PrintSelf(ostream &os, vtkIndent indent);

    vtkSetMacro(Tolerance, double);
    vtkGetMacro(Tolerance, double);
    vtkSetMacro(Scaling, double);
    vtkGetMacro(Scaling, double);

    vtkGetStringMacro(DownFileName);
    vtkSetStringMacro(DownFileName);

    vtkGetStringMacro(TopFileName);
    vtkSetStringMacro(TopFileName);

protected:

    vtkHhCalculator();
    ~vtkHhCalculator();

    virtual int RequestData(vtkInformation* request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);

    void loadPoints(const char* FileName, vtkSmartPointer<vtkPoints> sortedPoints);
    int UpdateLocators();

private:
    vtkHhCalculator(const vtkHhCalculator&);
    bool operator=(const vtkHhCalculator&);

    char* DownFileName;
    char* TopFileName;
    //
    char* LastLineAFile;
    char* LastLineBFile;
    double Tolerance;
    double Scaling;
    //
    vtkSmartPointer<vtkPolyData> lineA;
    vtkSmartPointer<vtkPolyData> lineB;
    //
    vtkSmartPointer<vtkDoubleArray> NormalsLineA;
    vtkSmartPointer<vtkDoubleArray> NormalsLineB;
};

#endif
