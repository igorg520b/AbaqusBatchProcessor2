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

    std::string rhitaD[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\nbd_1_1.msh)",
                         R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\nbd_2_1.msh)",
                         R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\nbd_3_1.msh)",
                         R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\nbd_4_1.msh)",
                         R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\nbd_5_1.msh)"};

    std::string rhitaL1[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl1_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl2_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl3_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl4_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl5_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl6_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl7_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl8_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl9_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL1\nbl10_1.msh)"};

    std::string rhitaL2[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl1_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl2_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl3_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl4_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl5_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl6_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl7_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl8_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl9_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaL2\nbl10_2.msh)"};



    std::string rhitaH1[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh1_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh2_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh3_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh4_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh5_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh6_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh7_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh8_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh9_1.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH1\nbh10_1.msh)"};

    std::string rhitaH2[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh1_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh2_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh3_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh4_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh5_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh6_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh7_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh8_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh9_2.msh)",
                           R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\rhitaH2\nbh10_2.msh)"};



    std::string keelL0[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn1.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn2.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn3.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn4.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn5.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn6.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn7.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn8.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn9.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL0\yn10.msh)"};

    std::string keelL1[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns1.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns2.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns3.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns4.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns5.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns6.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns7.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns8.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns9.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelL1\yns10.msh)"};

    std::string keelH0[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn1.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn2.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn3.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn4.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn5.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn6.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn7.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn8.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn9.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH0\wn10.msh)"};

    std::string keelH1[] {R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns1.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns2.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns3.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns4.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns5.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns6.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns7.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns8.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns9.msh)",
                          R"(C:\Users\s\Projects\GitHub\AbaqusBatchProcessor2\meshes\keelH1\wns10.msh)"};

    if(!QDir("results").exists()) QDir().mkdir("results");

    std::ofstream sjg, sje,sjexp;
    sjg.open("results\\jg.bat", std::ios_base::trunc|std::ios_base::out);
    sje.open("results\\jxc.bat", std::ios_base::trunc|std::ios_base::out);
    sjexp.open("results\\jexp.bat", std::ios_base::trunc|std::ios_base::out);

/*
    for(std::string s : rhitaD)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 1.0, 1.060325, 9);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }
    */


    // RHITA block
    for(std::string s : rhitaL1)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 1.0, 1.060325, 1);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(std::string s : rhitaL2)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 2.0, 1.060325, 2);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(std::string s : rhitaH1)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 1.0, 2.060325, 3);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(std::string s : rhitaH2)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 2.0, 2.060325, 4);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    // Keel
    for(std::string s : keelL0)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 5.0, 0, 5);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(std::string s : keelL1)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 6.0, 0, 6);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(std::string s : keelH0)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 5.0, 0, 7);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    for(std::string s : keelH1)
    {
        Generator g;
        g.LoadFromFileWithCrop(s, 6.0, 0, 8);
        g.AddEntryToJobGeneratorBat(sjg);
        g.AddEntryToJobExecutionBat(sje);
        g.AddEntryToBinaryExportBat(sjexp);
    }

    sjg.close();
    sje.close();
    sjexp.close();
}
