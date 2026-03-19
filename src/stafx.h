#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include <chrono>

#include <BSplCLib.hxx> 
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepLib.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeEdge2d.hxx>

#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomConvert.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomLib_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>


#include <gp_Pnt.hxx>
#include <gp_Circ.hxx>

#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <IGESControl_Reader.hxx>
#include <Interface_Static.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IntPatch_Polyhedron.hxx>

#include <math_Jacobi.hxx>

#include <TColgp_HArray1OfPnt.hxx>
#include <TColStd_HSequenceOfTransient.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_HSequenceOfShape.hxx>

#include <STEPControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <ShapeFix_Wireframe.hxx>

#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepProj_Projection.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>

#include <BRepGProp_Face.hxx>
#include <CPnts_AbscissaPoint.hxx>

#include <ShapeAnalysis_Edge.hxx>
#include <ShapeAnalysis_WireOrder.hxx>

#include <Eigen/Dense>

#include <numeric>
