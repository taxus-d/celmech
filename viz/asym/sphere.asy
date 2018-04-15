import three;
import solids;
import graph3;
unitsize(2cm);

settings.outformat="png";
settings.render=16;
settings.prc=false;

string imgdir="../../img/asym/";
string image_name = "sphere-small-transp";
string filename=imgdir + image_name;

currentprojection=perspective(
    camera=(4.2,-1.0,2.7),
    up=(0,0,1),
    target=(1,2.5,-2),
//     target=(-0.75768264256308,-0.164823812203309,0.16440107557311),
    angle=30//43.8077737959766
);
// currentlight.background = black + opacity(0.0);
real RE=1, RS=0.7, inc=100, lat=45, lon=90, tlat=50, tlon=100;

material spheremat = material(
    diffusepen=gray(0.5)+opacity(0.75),
    emissivepen=gray(0.3)
);
pen spheremesh = gray(0.35) + linewidth(1.3pt);
revolution shapesp=sphere((0,0,0), RE);
//Added specification for mesh pen
draw(surface(shapesp), surfacepen=spheremat, meshpen=spheremesh);

file in=input("../../data/orbit-ideal.dat").line().word();
real[][] a = in.dimension(1000,1000,1000);
real[][] data = transpose(a);


// real[] t=data[0];
real[] x=data[1];
real[] y=data[2];
real[] z=data[3];
draw(graph(x,y,z), red+linewidth(5pt));


pen collpen = blue+linewidth(3pt);
pen collpointpen = blue+linewidth(10pt);
triple[] collisions={(-1,0,0),(1/2,sqrt(3)/2,0),(1/2, -sqrt(3)/2,0)};
arrowbar3 collarrow = Arrow3(emissive(blue),size=20);

 //HookHead2, emissive=blue);
draw((0,0,0) -- 3*collisions[0], collpen, collarrow);
draw((0,0,0) -- 2*collisions[1], collpen,collarrow );
draw((0,0,0) -- 2*collisions[2], collpen,collarrow);


// dot((0,0,0), collpointpen);
dot(collisions[0], collpointpen);
dot(collisions[1], collpointpen);
dot(collisions[2], collpointpen);


shipout(filename);
// vim:ft=asy
