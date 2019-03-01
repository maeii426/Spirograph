// Global resolution
$fs = 0.5;  // Don't generate smaller facets than 0.1 mm
$fa = 5;    // Don't generate larger angles than 5 degrees

sharpie = 5;
radiusOuter = 9;
height = 8;

penHolder();

module penHolder() {
    
    union() {
    union() {
        difference(){
            difference() {
                cylinder(height, radiusOuter, radiusOuter, true);
                translate([0,0,-0.5]) cylinder(height+2, sharpie, sharpie, true);
            }
            cutOut(radiusOuter/2,height);
        }
        translate([0,0,0]) clasps(radiusOuter/2,height);
    }
    translate([-radiusOuter*2.25,0,0]) handle(radiusOuter,height);
    }
}

//Should be large enough to allow some standard sized
//screw through- use two washers to tighten around pen
module clasps(radius,height) {
    union() {
       translate([0,2*radius/5,0]) sideWall(radius,height);
       translate([0,-6*radius/5,0]) sideWall(radius,height);
    }
}

module handle(radius,height) {
CutPoints = [
  [  0,  -3*radius/5,  -height/2 ],  //0
  [ radius*1.5,  -3*radius/5,  -height/2 ],  //1
  [ radius*1.5,  3*radius/5,  -height/2 ],  //2
  [  0,  3*radius/5,  -height/2 ],  //3
  [  0,  -3*radius/5,  height/2 ],  //4
  [ radius*1.5,  -3*radius/5,  height/2 ],  //5
  [ radius*1.5,  3*radius/5,  height/2 ],  //6
  [  0,  3*radius/5,  height/2 ]]; //7
  
CutFaces = [
  [0,1,2,3],  // bottom
  [4,5,1,0],  // front
  [7,6,5,4],  // top
  [5,6,2,1],  // right
  [6,7,3,2],  // back
  [7,4,0,3]]; // left

difference() { 
    polyhedron(CutPoints, CutFaces);
    translate([radius*0.5,0,0]) cylinder(height*2, 2.5, 2.5, true);
}
}

module sideWall(radius,height) {
WallPoints = [
  [  sharpie,  0,  -height/2 ],  //0
  [ radius*3+3,  0,  -height/2 ],  //1
  [ radius*3+3,  radius*4/5,  -height/2 ],  //2
  [  sharpie+3,  radius*4/5,  -height/2 ],  //3
  [  sharpie,  0,  height/2 ],  //4
  [ radius*3+3,  0,  height/2 ],  //5
  [ radius*3+3,  radius*4/5,  height/2 ],  //6
  [  sharpie,  radius*4/5,  height/2 ]]; //7
  
WallFaces = [
  [0,1,2,3],  // bottom
  [4,5,1,0],  // front
  [7,6,5,4],  // top
  [5,6,2,1],  // right
  [6,7,3,2],  // back
  [7,4,0,3]]; // left
  
difference() {  
    polyhedron(WallPoints, WallFaces);
    rotate([90,0,0]) translate([radiusOuter+3,0,-2*radius/5]) cylinder(height*2, 2.5, 2.5, true);
}
}

module cutOut(radius,height) {
CutPoints = [
  [  0,  -2*radius/5,  -height ],  //0
  [ radius*2,  -2*radius/5,  -height ],  //1
  [ radius*2,  2*radius/5,  -height ],  //2
  [  0,  2*radius/5,  -height ],  //3
  [  0,  -2*radius/5,  height ],  //4
  [ radius*2,  -2*radius/5,  height ],  //5
  [ radius*2,  2*radius/5,  height ],  //6
  [  0,  2*radius/5,  height ]]; //7
  
CutFaces = [
  [0,1,2,3],  // bottom
  [4,5,1,0],  // front
  [7,6,5,4],  // top
  [5,6,2,1],  // right
  [6,7,3,2],  // back
  [7,4,0,3]]; // left
    
polyhedron(CutPoints, CutFaces);
}