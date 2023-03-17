#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDebug>
#include <QFileInfo>
#include <QDir>

#include "gmsh.h"
#include "generator.h"

#include <string>
#include <fstream>


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


    if(!QDir("results").exists()) QDir().mkdir("results");
    if(!QDir("inp").exists()) QDir().mkdir("inp");
    if(!QDir("exportscripts").exists()) QDir().mkdir("exportscripts");

    std::ofstream sjg, sje,sjexp;
    sjg.open("results\\jg.bat", std::ios_base::trunc|std::ios_base::out);
    sje.open("inp\\jxc.bat", std::ios_base::trunc|std::ios_base::out);
    sjexp.open("results\\jexp.bat", std::ios_base::trunc|std::ios_base::out);

    // RHITA block
    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(rhitaL1[i], 1.0, 1.060325, 1);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(rhitaL2[i], 2.0, 1.060325, 2);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(rhitaH1[i], 1.0, 2.060325, 3);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }


    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(rhitaH2[i], 2.0, 2.060325, 4);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    // Keel
    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(keelL0[i], 5.0, 0, 5);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(keelL1[i], 6.0, 0, 6);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(keelH0[i], 5.0, 0, 7);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(int i=0;i<5;i++)
    {
        Generator g;
        g.LoadFromFileWithCrop(keelH1[i], 6.0, 0, 8);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    sjg.close();
    sje.close();
    sjexp.close();
}
