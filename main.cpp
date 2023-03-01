#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDebug>
#include <QFileInfo>

#include "gmsh.h"
#include "generator.h"

#include <string>


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QCoreApplication::setApplicationName("RHITA geometry generator");
    QCoreApplication::setApplicationVersion("1.0");

    QCommandLineParser parser;
    parser.setApplicationDescription("Batch generator of Abaqus .inp files");
    parser.addHelpOption();
    parser.addVersionOption();
    //    parser.addPositionalArgument("source", QCoreApplication::translate("main", "Configuration file."));

    std::string rhitaL1[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL1\\nbl1_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL1\\nbl2_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL1\\nbl3_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL1\\nbl4_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL1\\nbl5_1.msh"};

    std::string rhitaL2[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL2\\nbl1_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL2\\nbl2_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL2\\nbl3_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL2\\nbl4_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL2\\nbl5_2.msh"};

    std::string rhitaH1[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH1\\nbh1_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH1\\nbh2_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH1\\nbh3_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH1\\nbh4_1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH1\\nbh5_1.msh"};

    std::string rhitaH2[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH2\\nbh1_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH2\\nbh2_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH2\\nbh3_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH2\\nbh4_2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaH2\\nbh5_2.msh"};


    std::string keelL0[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL0\\yn1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL0\\yn2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL0\\yn3.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL0\\yn4.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL0\\yn5.msh"};

    std::string keelL1[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL1\\yns1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL1\\yns2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL1\\yns3.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL1\\yns4.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelL1\\yns5.msh"};

    std::string keelH0[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH0\\wn1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH0\\wn2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH0\\wn3.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH0\\wn4.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH0\\wn5.msh"};

    std::string keelH1[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH1\\wns1.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH1\\wns2.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH1\\wns3.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH1\\wns4.msh",
                          "C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\keelH1\\wns5.msh"};


    for(int i=0;i<5;i++)
    {
        Generator g;
//        g.LoadFromFileWithCrop(rhitaL1[i], 1.0, 1.111925, 1);
//        g.LoadFromFileWithCrop(rhitaL2[i], 2.0, 1.111925, 2);
//        g.LoadFromFileWithCrop(rhitaH1[i], 1.0, 2.111925, 3);
//        g.LoadFromFileWithCrop(rhitaH2[i], 2.0, 2.111925, 4);
//        g.LoadFromFileWithCrop(keelL0[i], 5.0, 0, 5);
//        g.LoadFromFileWithCrop(keelL1[i], 6.0, 0, 6);
//        g.LoadFromFileWithCrop(keelH0[i], 5.0, 0, 7);
        g.LoadFromFileWithCrop(keelH1[i], 6.0, 0, 8);
    }



}
