#include "CurveAndSurfaceTestCases.h"
#include "AlgAutoCurFPIA.h"
#include "AlgAutoSurfFPIA.h"

void CurveAndSurfaceTestCases::surfaceCase3_Car()
{
    Handle(Geom_BSplineSurface) surface;
    readIGESFileToSurface("111-surface.igs", surface);
    AlgSurfFPIAByAdjustingControlPoint alg2;
    alg2.Init(surface);
    int poleNum = surface->NbUPoles() * surface->NbVPoles();
    std::vector<size_t> indexList(poleNum);
    std::vector<double> weights(poleNum, 0.1);
    for (int i = 0; i < poleNum; i++)
    {
        indexList[i] = i;
        if (i == 12 || i == 13 || i == 14 || i== 19) weights[i] = 0.18;
    }
    alg2.setWeights(indexList, weights);
    alg2.execute(1000, 1e-4);
    //SurfaceToIgs(alg2.getResult(), "..\\Surface-Face-pia.igs");
}

void CurveAndSurfaceTestCases::curveCase1_G()
{
    //读取数据拟合成曲线
    Handle(Geom_BSplineCurve) curve;
    readFileToCurve("..\\test data\\G-shape537.cur", curve);
   
    //PIA
    auto algTest = std::make_shared<AlgCurFPIAByAdjustingControlPoint>();
    algTest->Init(curve, 2);
    //设置具体的光顺权重
    int poleNum = curve->NbPoles();
    std::vector<size_t> indexList(poleNum);
    std::vector<double> weights(poleNum,1e-6);
    for (int i = 0; i < poleNum; i++)
    {
        indexList[i] = i;
        if (i > 0 && i < 6)
            weights[i] = 1e-7;
        if (i > 34 && i < 40)
            weights[i] = 2e-6;
        if (i > 62 && i < 68)
            weights[i] = 2e-6;
    }
    algTest->execute(indexList, weights, 800);
    //algTest->execute(1000, 1e-4);
    //CurveToIgs(algTest->getResult(), "..\\test data\\G-res.igs");
}

void CurveAndSurfaceTestCases::curveCase2_curve()
{
    Handle(Geom_BSplineCurve) curve;
    readFileToCurve("..\\test data\\dolphin.cur", curve);

    auto algTest = std::make_shared<AlgCurFPIAByAdjustingControlPoint>();
    algTest->Init(curve, 2);
    //设置具体的光顺权重
    int poleNum = curve->NbPoles();
    std::vector<size_t> indexList(poleNum);
    std::vector<double> weights(poleNum, 1e-6);
    for (int i = 0; i < poleNum; i++)
    {
        indexList[i] = i;
        if (i > 24 && i < 33)
            weights[i] = 1e-5;
        if (i > 45 && i < 58)
            weights[i] = 8e-6;
    }
    algTest->execute(indexList, weights, 1000);
}

void CurveAndSurfaceTestCases::curveCase3()
{
    int num = 100;
    double a = 2.0;
    double b = 1.5;
    std::vector<gp_Pnt> data(num);
    for (int i = 0; i < num; i++)
    {
        double t = i * 0.05;
        double r = a + b * t;
        double x = r * cos(t);
        double y = r * sin(t);
        gp_Pnt point{ x,y,0 };
        data[i] = point;
        //std::cout << x << " " << y << std::endl;
    }
    //std::cout << "==================================\n";
    auto newData = addGaussianNoiseToPoint(data, 30, 60, 0, 0.12);
    newData = addGaussianNoiseToPoint(newData, 70, 80, 0, 0.05);

    AlgPointFitCurve alg1;
    alg1.Init(newData, 30, 1);
    Handle(Geom_BSplineCurve)  theCurve1 = alg1.getResult();
    
    AlgPointFitCurve alg2;
    alg2.Init(newData, 30, 2);
    Handle(Geom_BSplineCurve)  theCurve2 = alg2.getResult();

    AlgPointFitCurve alg3;
    alg3.Init(newData, 30, 3);
    Handle(Geom_BSplineCurve)  theCurve3 = alg3.getResult();
}

void CurveAndSurfaceTestCases::curveAutoCase()
{
    Handle(Geom_BSplineCurve) curve;
    readFileToCurve("..\\test data\\cat.cur", curve);
    AlgAutoCurFPIA alg(curve, 2);

    std::vector<double> weights = {1e-6,1e-6,1e-6};
    alg.Fairing(3, weights);
    curve = alg.getResult();

    //std::vector<double> weights = {5e-5, 5e-5,5e-5,5e-5,1e-6,1e-6,1e-6,1e-6 };
    //alg.Fairing(8, weights);
    //curve = alg.getResult();
}

void CurveAndSurfaceTestCases::surfaceAutoCase()
{
    Handle(Geom_BSplineSurface) surface;
    readIGESFileToSurface("..\\555-surface.igs", surface);

    AlgAutoSurfFPIA alg(surface, 2);
    std::vector<double> weights = { 0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005 };
    alg.Fairing(13, weights);
    surface = alg.getResult();
}

void CurveAndSurfaceTestCases::SurfaceToIgs(const Handle(Geom_BSplineSurface)& surface, Standard_CString theFileName)
{
    BRepBuilderAPI_MakeFace anFaceMaker(surface, 1e-3);
    TopoDS_Shape aShape = anFaceMaker.Face();

    IGESControl_Writer aWriter;
    Standard_Boolean status = Standard_True;

    status = status && aWriter.AddShape(aShape);
    if (status == Standard_False)
    {
        return;
    }
    status = aWriter.Write(theFileName);
}

void CurveAndSurfaceTestCases::CurveToIgs(const Handle(Geom_BSplineCurve)& fairingCurve, Standard_CString theFileName)
{
    BRepBuilderAPI_MakeEdge anEdgeMaker(fairingCurve);
    TopoDS_Shape aShape = anEdgeMaker.Edge();

    IGESControl_Writer aWriter;
    Standard_Boolean status = Standard_True;

    status = status && aWriter.AddShape(aShape);
    if (status == Standard_False)
    {
        return;
    }
    status = aWriter.Write(theFileName);
}

void CurveAndSurfaceTestCases::PointsToPLY(const std::vector<gp_Pnt>& points, const std::string& filename)
{
    std::ofstream file(filename);
    file << "ply\nformat ascii 1.0\n";
    file << "element vertex " << points.size() << "\n";
    file << "property float x\nproperty float y\nproperty float z\n";
    file << "end_header\n";

    for (const auto& p : points) {
        file << p.X() << " " << p.Y() << " " << p.Z() << "\n";
    }
    std::cout << "Saved " << points.size() << " points to " << filename << std::endl;
}

bool CurveAndSurfaceTestCases::readFileToData(const Standard_CString& theFileName, std::vector<std::vector<gp_Pnt>>& theData)
{
    std::ifstream readFile;
    readFile.open(theFileName, std::ios::in);

    char s[50];
    std::string str;

    int row = 0, col = 0;

    std::random_device rd;  // 真随机数生成器
    std::mt19937 gen(rd()); // Mersenne Twister 引擎
    std::uniform_int_distribution<> dis(0, 99);  // 生成 0 到 99 之间的整数

    if (readFile.is_open())
    {
        readFile >> str; // # 
        readFile >> str; // REM
        readFile >> str; // format
        readFile >> str; // file

        readFile >> str;
        row = stoi(str);

        readFile >> str;
        col = stoi(str);

        int datanum = row * col;
        std::cout << "数据规模：" << row << " " << col << std::endl;
        theData.resize(row, std::vector<gp_Pnt>(col));

        int p = 0;
        while (p < datanum)
        {
            readFile >> str; // v

            readFile >> str; // x
            str = str.substr(0, 10);
            double x = stod(str);

            readFile >> str; // y
            str = str.substr(0, 10);
            double y = stod(str);

            readFile >> str; // z
            str = str.substr(0, 10);
            double z = stod(str);

            int i1 = p / col;
            int i2 = p % col;
            theData[i1][i2] = gp_Pnt(x, y, z);
            p++;
        }
    }
    else
    {
        std::cout << "文件打开失败\n";
    }

    readFile.close();
}

void CurveAndSurfaceTestCases::readFileToCurve(std::string fileName, Handle(Geom_BSplineCurve)& theCurve)
{
    std::ifstream readFile;
    readFile.open(fileName, std::ios::in);

    char s[50];
    std::string str;

    int datanum = 0;
    int dim = 2;

    std::random_device rd;  // 真随机数生成器
    std::mt19937 gen(rd()); // Mersenne Twister 引擎
    std::uniform_int_distribution<> dis(0, 99);  // 生成 0 到 99 之间的整数

    std::vector<gp_Pnt> theData;

    if (readFile.is_open())
    {
        readFile >> str;
        readFile >> str;
        datanum = stoi(str);
        std::cout << datanum << std::endl;
        int p = 0;
        theData.resize(datanum);
        //int index = 0;
        while (p < datanum)
        {

            readFile >> str;//x
            str = str.substr(0, 10);
            double x = stod(str);

            readFile >> str;//y
            str = str.substr(0, 10);
            double y = stod(str);

            theData[p] = gp_Pnt(x, y, 0);
            p++;
        }
    }
    else
    {
        std::cout << "文件打开失败\n";
    }

    AlgPointFitCurve alg;
    alg.Init(theData, 100, 3);
    theCurve = alg.getResult();

    readFile.close();
}

bool CurveAndSurfaceTestCases::readIGESFileToCurve(const Standard_CString& theFileName, Handle(Geom_BSplineCurve)& theCurve)
{
    IGESControl_Reader aReader;
    IFSelect_ReturnStatus aStatus = aReader.ReadFile(theFileName);
    if (aStatus != IFSelect_RetDone)
    {
        return false;
    }

    aReader.TransferRoots();
    Standard_Integer aNumbs = aReader.NbShapes();

    //std::cout << aNumbs << std::endl;
    if (aNumbs == 1)
    {

        TopoDS_Shape aShape = aReader.Shape(1);
        TopExp_Explorer aEdgeExp(aShape, TopAbs_EDGE);

        TopoDS_Edge aTempEdge = TopoDS::Edge(aEdgeExp.Current());
        BRepAdaptor_Curve aTempBRepAdaptor(aTempEdge);

        GeomAdaptor_Curve aTempGeomAdaptor = aTempBRepAdaptor.Curve();
        Handle(Geom_Curve) curve = aTempGeomAdaptor.Curve();
        Handle(Geom_BSplineCurve) bcurve = GeomConvert::CurveToBSplineCurve(curve);
        //std::cout << bcurve->Degree() << std::endl;
        //std::cout << bcurve->NbPoles() << std::endl;
        theCurve = bcurve;
        //checkCurve(theCurve);
    }
    return true;
}

bool CurveAndSurfaceTestCases::readIGESFileToSurface(const Standard_CString& theFileName, Handle(Geom_BSplineSurface)& theSurface)
{
    IGESControl_Reader aReader;
    IFSelect_ReturnStatus aStatus = aReader.ReadFile(theFileName);
    if (aStatus != IFSelect_RetDone)
    {
        return false;
    }

    aReader.TransferRoots();
    Standard_Integer aNumbs = aReader.NbShapes();

    std::cout << aNumbs << std::endl;
    if (aNumbs == 1)
    {

        TopExp_Explorer anExplorer(aReader.Shape(), TopAbs_FACE);
        const TopoDS_Face& aFace = TopoDS::Face(anExplorer.Current());
        BRepAdaptor_Surface aBRepAdaptor(aFace, Standard_True);
        GeomAdaptor_Surface aGeomAdaptor = aBRepAdaptor.Surface();
        if (aGeomAdaptor.GetType() == GeomAbs_BSplineSurface)
        {
            theSurface = GeomConvert::SurfaceToBSplineSurface(aGeomAdaptor.Surface());
        }
    }
    return true;
}

void CurveAndSurfaceTestCases::createSpiralCurve(Handle(Geom_BSplineCurve)& theCurve)
{
    int num = 100;
    double a = 2.0;
    double b = 1.5;
    std::vector<gp_Pnt> data(num);
    for (int i = 0; i < num; i++)
    {
        double t = i * 0.05;
        double r = a + b * t;
        double x = r * cos(t);
        double y = r * sin(t);
        gp_Pnt point{ x,y,0 };
        data[i] = point;
        //std::cout << x << " " << y << std::endl;
    }
    //std::cout << "==================================\n";
    auto newData = addGaussianNoiseToPoint(data, 30, 60, 0, 0.12);
    newData = addGaussianNoiseToPoint(newData, 70, 80, 0, 0.05);

    AlgPointFitCurve alg;
    alg.Init(newData, 30, 3);
    theCurve = alg.getResult();
    
}

std::vector<gp_Pnt> CurveAndSurfaceTestCases::addGaussianNoiseToPoint(const std::vector<gp_Pnt>& data, int begin, int end, double mean, double stddev)
{
    std::vector<gp_Pnt> noisyData(data);

    // 创建随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(mean, stddev);

    // 添加噪声
    /*for (gp_Pnt& p : noisyData) {
        double x = p.X() + d(gen);
        double y = p.Y() + d(gen);
        p.SetX(x);
        p.SetY(y);
    }*/

    for (int i = begin; i <= end; i++)
    {
        gp_Pnt p = noisyData[i];
        double x = p.X() + d(gen);
        double y = p.Y() + d(gen);
        p.SetX(x);
        p.SetY(y);
        noisyData[i] = p;
    }

    return noisyData;
}

std::vector<std::vector<gp_Pnt>> CurveAndSurfaceTestCases::addGaussianNoiseToSurfacePoint(const std::vector<std::vector<gp_Pnt>>& data, int begin1, int end1, int begin2, int end2, double mean, double stddev)
{
    std::vector<std::vector<gp_Pnt>> noisyData(data);

    // 创建随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(mean, stddev);

    // 添加噪声
    for (int i = begin1; i <= end1; i++)
    {
        for (int j = begin2; j <= end2; j++)
        {
            gp_Pnt p = noisyData[i][j];
            double x = p.X() + d(gen);
            double y = p.Y() + d(gen);
            double z = p.Z() + d(gen);

            p.SetX(x);
            p.SetY(y);
            p.SetZ(z);

            noisyData[i][j] = p;
        }

    }

    return noisyData;
}

void CurveAndSurfaceTestCases::checkCurve(const Handle(Geom_BSplineCurve)& theCurve)
{
    if (theCurve.IsNull())
    {
        std::cout << "曲线为空\n";
        return;
    }

    std::cout << "曲线次数" << theCurve->Degree() << std::endl;
    std::cout << "控制点个数" << theCurve->NbPoles() << std::endl;

    auto poles = theCurve->Poles();
    for (auto& p : poles)
    {
        std::cout << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")\n";
    }
}

void CurveAndSurfaceTestCases::checkSurface(const Handle(Geom_BSplineSurface)& theSurface)
{
    if (theSurface.IsNull())
    {
        std::cout << "曲面为空\n";
        return;
    }

    std::cout << "曲面次数" << theSurface->UDegree() << " " << theSurface->VDegree() << std::endl;
    std::cout << "控制点个数" << theSurface->NbUPoles() << " " << theSurface->NbVPoles() << std::endl;

    auto poles = theSurface->Poles();
    for (auto& p : poles)
    {
        std::cout << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")\n";
    }
}
