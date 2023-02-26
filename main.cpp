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

    std::string files[] {"C:\\Users\\s\\Projects\\GitHub\\AbaqusBatchProcessor2\\meshes\\rhitaL1\\nbl1_1.msh"};

    Generator g;
    g.LoadFromFileWithCrop(files[0], 1.0, 1.111925, 1);



}
