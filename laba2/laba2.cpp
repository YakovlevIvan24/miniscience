#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;
double tau = 0.01;
// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Некая величина, в попугаях
    double temp;
    // Скорость
    double vx;
    double vy;
    double vz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), temp(25), vx(0.0), vy(0.0), vz(0.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double temp, double vx, double vy, double vz) 
            : x(x), y(y), z(z), temp(temp), vx(vx), vy(vy), vz(vz)
    {
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double t, unsigned int flag,double top) {
        double dt1 = 0;
        double dt2 = 0;
        switch (flag)
        {

            case 0:
                vy = 1;
            break;
            case 1:
                if(y < top -2.665){ vy = 1;}
                else {vy = 3;};
            break;
            case 2:
                dt1 = t - 20;
                dt2 = t -80;
//                if(y < top -1.095 & y >= top -2.665 ){vy =2;}
//                else if(y < top -2.665){ vy = 1;}
//                else{vy = 5;}
                if(y > top -1.04){vy = 5;}
                else if (y < top -1.095 - 2*dt2/100 & y > top -2.67 - 2*dt2/100 ){vy =3;}
                else if(y <= top -2.67-2*dt2/100-2*(dt2+dt1)/100){ vy = 1;}

            break;
        }
        vx = 0;
        vz = 0;
        x += vx*tau;
        y += vy*tau;
        z += vz*tau;
        double lambda = 3;
        double v = 0.03;
        double nu = v/lambda;
        double omega = 2*M_PI*nu;
        double k = 2*M_PI/lambda;
        double rho = sqrt(pow((x),2) + pow(top-y,2) + pow(z,2));
        temp =  100000000000000/rho*cos(omega*t - k*rho) + 100000000000000/rho*cos(omega*t + k*rho);
    }

};

// Класс элемента сетки
class Element
{
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {
        double top = 0;
        for(int i = 0; i< nodes.size()/3; i++)
        {
            top = max(top, nodesCoords[i*3 + 1]);
        }

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            // Модельная скалярная величина распределена как-то вот так
            double temp = 20*pow(pointY, 2)/(pow(top,2));
            nodes[i] = CalcNode(pointX, pointY, pointZ, temp, 0.0, 0.0, 0.0);
        }

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double t, unsigned int flag, double top) {
        
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            nodes[i].move(t,flag,top);
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("temp");
        auto velY = vtkSmartPointer<vtkDoubleArray>::New();
        velY->SetName("Y velocity");
        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
            vel->InsertNextTuple(_vel);
            velY->InsertNextValue(nodes[i].vy);
            // И значение скалярного поля тоже
            smth->InsertNextValue(nodes[i].temp);
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);
        unstructuredGrid->GetPointData()->AddArray(velY);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "ROCKET-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
    double getTop(){
        double top = 0;
        for(int i = 0; i< nodes.size()/3; i++)
        {
            top = max(top, nodes[i*3 + 1].y);
        }
        return top;
    }
};

int main()
{
    // Шаг точек по пространству
    double h = 4.0;
    // Шаг по времени

    const unsigned int GMSH_TETR_CODE = 4;

    // Теперь придётся немного упороться:
    // (а) построением сетки средствами gmsh,
    // (б) извлечением данных этой сетки в свой код.
    gmsh::initialize();

//   gmsh::model::add("ROCKET1");

//   // Let's merge an STL mesh that we would like to remesh (from the parent
//   // directory):
//   try {
//     gmsh::merge("ROCKET4.stl");
//   } catch(...) {
//     gmsh::logger::write("Could not load STL mesh: bye!");
//     gmsh::finalize();
//     return 0;
//   }
//     gmsh::option::setNumber("Mesh.Algorithm3D",1);
//   // We first classify ("color") the surfaces by splitting the original surface
//   // along sharp geometrical features. This will create new discrete surfaces,
//   // curves and points.

//     gmsh::model::mesh::classifySurfaces(0 * M_PI / 180., true, true);

//     gmsh::model::mesh::createGeometry();
//   // Angle between two triangles above which an edge is considered as sharp:


//   // For complex geometries, patches can be too complex, too elongated or too
//   // large to be parametrized; setting the following option will force the
//   // creation of patches that are amenable to reparametrization

//   // Create a geometry for all the discrete curves and surfaces in the mesh, b

//   // Create a volume from all the surfaces
//   std::vector<std::pair<int, int> > s;
//   gmsh::model::getEntities(s, 2);
//   std::vector<int> sl;
//   for(auto surf : s) sl.push_back(surf.second);
//   int l = gmsh::model::geo::addSurfaceLoop(sl);
//   gmsh::model::geo::addVolume({l});

//   gmsh::model::geo::synchronize();

//   // We specify element sizes imposed by a size field, just because we can :-)
//   gmsh::model::mesh::generate(3);
    gmsh::open("FINAL_ROCKET.msh");
    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    // TODO: неплохо бы полноценно данные сетки проверять, да

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);
    double top = mesh.getTop();
    gmsh::finalize();
    unsigned int flag = 0;
    for(int i = 0; i < 240;i++)
    {
        if(i == 20) {flag = 1;}
        else if(i == 80) {flag = 2;};
        top = mesh.getTop();
        mesh.doTimeStep(i,flag,top);
        mesh.snapshot(i);
    }
    return 0;
}
