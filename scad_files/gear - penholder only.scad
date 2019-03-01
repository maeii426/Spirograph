//////////////////////////////////////////////////////////////////////////////////////////////
// Public Domain Parametric Involute Spur Gear (and involute helical gear and involute rack)
// version 1.1
// by Leemon Baird, 2011, Leemon@Leemon.com
//http://www.thingiverse.com/thing:5505
//
// This file is public domain.  Use it for any purpose, including commercial
// applications.  Attribution would be nice, but is not required.  There is
// no warranty of any kind, including its correctness, usefulness, or safety.
// 
// This is parameterized involute spur (or helical) gear.  It is much simpler and less powerful than
// others on Thingiverse.  But it is public domain.  I implemented it from scratch from the 
// descriptions and equations on Wikipedia and the web, using Mathematica for calculations and testing,
// and I now release it into the public domain.
//
//		http://en.wikipedia.org/wiki/Involute_gear
//		http://en.wikipedia.org/wiki/Gear
//		http://en.wikipedia.org/wiki/List_of_gear_nomenclature
//		http://gtrebaol.free.fr/doc/catia/spur_gear.html
//		http://www.cs.cmu.edu/~rapidproto/mechanisms/chpt7.html
//
// The module gear() gives an involute spur gear, with reasonable defaults for all the parameters.
// Normally, you should just choose the first 4 parameters, and let the rest be default values.
// The module gear() gives a gear in the XY plane, centered on the origin, with one tooth centered on
// the positive Y axis.  The various functions below it take the same parameters, and return various
// measurements for the gear.  The most important is pitch_radius, which tells how far apart to space
// gears that are meshing, and adendum_radius, which gives the size of the region filled by the gear.
// A gear has a "pitch circle", which is an invisible circle that cuts through the middle of each
// tooth (though not the exact center). In order for two gears to mesh, their pitch circles should 
// just touch.  So the distance between their centers should be pitch_radius() for one, plus pitch_radius() 
// for the other, which gives the radii of their pitch circles.
//
// In order for two gears to mesh, they must have the same mm_per_tooth and pressure_angle parameters.  
// mm_per_tooth gives the number of millimeters of arc around the pitch circle covered by one tooth and one
// space between teeth.  The pitch angle controls how flat or bulged the sides of the teeth are.  Common
// values include 14.5 degrees and 20 degrees, and occasionally 25.  Though I've seen 28 recommended for
// plastic gears. Larger numbers bulge out more, giving stronger teeth, so 28 degrees is the default here.
//
// The ratio of number_of_teeth for two meshing gears gives how many times one will make a full 
// revolution when the the other makes one full revolution.  If the two numbers are coprime (i.e. 
// are not both divisible by the same number greater than 1), then every tooth on one gear
// will meet every tooth on the other, for more even wear.  So coprime numbers of teeth are good.
//
// The module rack() gives a rack, which is a bar with teeth.  A rack can mesh with any
// gear that has the same mm_per_tooth and pressure_angle.
//
// Some terminology: 
// The outline of a gear is a smooth circle (the "pitch circle") which has mountains and valleys
// added so it is toothed.  So there is an inner circle (the "root circle") that touches the 
// base of all the teeth, an outer circle that touches the tips of all the teeth,
// and the invisible pitch circle in between them.  There is also a "base circle", which can be smaller than
// all three of the others, which controls the shape of the teeth.  The side of each tooth lies on the path 
// that the end of a string would follow if it were wrapped tightly around the base circle, then slowly unwound.  
// That shape is an "involute", which gives this type of gear its name.
//
//////////////////////////////////////////////////////////////////////////////////////////////

//An involute spur gear, with reasonable defaults for all the parameters.
//Normally, you should just choose the first 4 parameters, and let the rest be default values.
//Meshing gears must match in mm_per_tooth, pressure_angle, and twist,
//and be separated by the sum of their pitch radii, which can be found with pitch_radius().
module gear (
	mm_per_tooth    = 3,    //this is the "circular pitch", the circumference of the pitch circle divided by the number of teeth
	number_of_teeth = 11,   //total number of teeth around the entire perimeter
	thickness       = 6,    //thickness of gear in mm
	hole_diameter   = 3,    //diameter of the hole in the center, in mm
	twist           = 0,    //teeth rotate this many degrees from bottom of gear to top.  360 makes the gear a screw with each thread going around once
	teeth_to_hide   = 0,    //number of teeth to delete to make this only a fraction of a circle
	pressure_angle  = 28,   //Controls how straight or bulged the tooth sides are. In degrees.
	clearance       = 0.0,  //gap between top of a tooth on one gear and bottom of valley on a meshing gear (in millimeters)
	backlash        = 0.0   //gap between two meshing teeth, in the direction along the circumference of the pitch circle
) {
	assign(pi = 3.1415926)
	assign(p  = mm_per_tooth * number_of_teeth / pi / 2)  //radius of pitch circle
	assign(c  = p + mm_per_tooth / pi - clearance)        //radius of outer circle
	assign(b  = p*cos(pressure_angle))                    //radius of base circle
	assign(r  = p-(c-p)-clearance)                        //radius of root circle
	assign(t  = mm_per_tooth/2-backlash/2)                //tooth thickness at pitch circle
	assign(k  = -iang(b, p) - t/2/p/pi*180) {             //angle to where involute meets base circle on each side of tooth
		union() {
			for (i = [0:number_of_teeth-teeth_to_hide-1] )
				rotate([0,0,i*360/number_of_teeth])
					linear_extrude(height = thickness, center = true, convexity = 10, twist = twist)
						polygon(
							points=[
								[0, -hole_diameter/10],
								polar(r, -181/number_of_teeth),
								polar(r, r<b ? k : -180/number_of_teeth),
								q7(0/5,r,b,c,k, 1),q7(1/5,r,b,c,k, 1),q7(2/5,r,b,c,k, 1),q7(3/5,r,b,c,k, 1),q7(4/5,r,b,c,k, 1),q7(5/5,r,b,c,k, 1),
								q7(5/5,r,b,c,k,-1),q7(4/5,r,b,c,k,-1),q7(3/5,r,b,c,k,-1),q7(2/5,r,b,c,k,-1),q7(1/5,r,b,c,k,-1),q7(0/5,r,b,c,k,-1),
								polar(r, r<b ? -k : 180/number_of_teeth),
								polar(r, 181/number_of_teeth)
							],
 							paths=[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]]
						);
			translate([0,0,-3]) cylinder(h=2*thickness+1, r=hole_diameter/2, center=true, $fn=20);
		}
	}
};	
//these 4 functions are used by gear
function polar(r,theta)   = r*[sin(theta), cos(theta)];                            //convert polar to cartesian coordinates
function iang(r1,r2)      = sqrt((r2/r1)*(r2/r1) - 1)/3.1415926*180 - acos(r1/r2); //unwind a string this many degrees to go from radius r1 to radius r2
function q7(f,r,b,r2,t,s) = q6(b,s,t,(1-f)*max(b,r)+f*r2);                         //radius a fraction f up the curved side of the tooth 
function q6(b,s,t,d)      = polar(d,s*(iang(b,d)+t));                              //point at radius d on the involute curve



//These 5 functions let the user find the derived dimensions of the gear.
//A gear fits within a circle of radius outer_radius, and two gears should have
//their centers separated by the sum of their pictch_radius.
function circular_pitch  (mm_per_tooth=3) = mm_per_tooth;                     //tooth density expressed as "circular pitch" in millimeters
function diametral_pitch (mm_per_tooth=3) = 3.1415926 / mm_per_tooth;         //tooth density expressed as "diametral pitch" in teeth per millimeter
function module_value    (mm_per_tooth=3) = mm_per_tooth / pi;                //tooth density expressed as "module" or "modulus" in millimeters
function pitch_radius    (mm_per_tooth=3,number_of_teeth=11) = mm_per_tooth * number_of_teeth / 3.1415926 / 2;
function outer_radius    (mm_per_tooth=3,number_of_teeth=11,clearance=0.1)    //The gear fits entirely within a cylinder of this radius.
	= mm_per_tooth*(1+number_of_teeth/2)/3.1415926  - clearance;              

//////////////////////////////////////////////////////////////////////////////////////////////
//example gear train.  
//Try it with OpenSCAD View/Animate command with 20 steps and 24 FPS.
//The gears will continue to be rotated to mesh correctly if you change the number of teeth.

n1 = 13; //red gear number of teeth
echo("num",ceil(n1*(3.6)));
n2 = ceil(n1*(3.6)); //green gear
n3 = ceil(n1*(1.5)); //extra gear
//n3 = 5;  //blue gear
//n4 = 20; //orange gear
//n5 = 8;  //gray rack
mm_per_tooth = 8; //all meshing gears need the same mm_per_tooth (and the same pressure_angle)
thickness    = 5;
vert_shift = thickness/2;
big_gear_thickness = thickness*1.2;
big_gear_vert_shift = big_gear_thickness/2;
hole         = 5; // diameter in mm
height       = 12;

d1 =pitch_radius(mm_per_tooth,n1);
d12=pitch_radius(mm_per_tooth,n1) + pitch_radius(mm_per_tooth,n2);
//d13=pitch_radius(mm_per_tooth,n1) + pitch_radius(mm_per_tooth,n3);
//d14=pitch_radius(mm_per_tooth,n1) + pitch_radius(mm_per_tooth,n4);

o_r1=outer_radius(mm_per_tooth,n1);
o_r2=outer_radius(mm_per_tooth,n2);
p_r1=pitch_radius(mm_per_tooth,n1);
p_r2=pitch_radius(mm_per_tooth,n2);

pin_radius    = (p_r1-(hole/2))/4.5;
pin_height   = o_r1*1.75;
pin_shift    = (hole/2)+(p_r1*(3/7));
canvas_thickness = big_gear_thickness/4;
bar_width = pin_radius*4;
bar_length = o_r1*11.76;
bar_hole_length = bar_length*.9;
support_radius = hole;
support_height = pin_height/2-thickness/2-1.5;

sharpie = 5;
radiusOuter = 9;
height = 8;
penHolder();


module penHolder() {
    union() {
        union() {
            union() {
                difference(){
                    difference() {
                        cylinder(height, radiusOuter, radiusOuter, true, $fn=40);
                        translate([0,0,-0.5]) cylinder(height+2, sharpie, sharpie, true, $fn=40);
                    }
                    cutOut(radiusOuter/2,height);
                }
                translate([0,0,0]) clasps(radiusOuter/2,height);
            }
        translate([-radiusOuter*2.25,0,0]) handle(radiusOuter,height);
        }
        // support under penholder
        difference () {
            translate([-thickness*1.69,0,-thickness]) cube(size=[thickness/2,thickness*2,thickness*1.5], center=true);
            translate([-thickness*1.4,0,-thickness*1.27]) cube(size=[bar_width,1,thickness*.75], center=true);
        }
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
  [ -radius*.25,  -3*radius/5,  -height/2 ],  //0
  [ radius*1.5,  -3*radius/5,  -height/2 ],  //1
  [ radius*1.5,  3*radius/5,  -height/2 ],  //2
  [ -radius*.25,  3*radius/5,  -height/2 ],  //3
  [ -radius*.25,  -3*radius/5,  height/2 ],  //4
  [ radius*1.5,  -3*radius/5,  height/2 ],  //5
  [ radius*1.5,  3*radius/5,  height/2 ],  //6
  [ -radius*.25,  3*radius/5,  height/2 ]]; //7
  
CutFaces = [
  [0,1,2,3],  // bottom
  [4,5,1,0],  // front
  [7,6,5,4],  // top
  [5,6,2,1],  // right
  [6,7,3,2],  // back
  [7,4,0,3]]; // left

difference() { 
    polyhedron(CutPoints, CutFaces);
    translate([radius*0.5,0,0]) cylinder(height*2, 2.5, 2.5, true, $fn=20);
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
    rotate([90,0,0]) translate([radiusOuter+3,0,-2*radius/5]) cylinder(height*2, 2.5, 2.5, true, $fn=20);
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

    
echo("outer radius small: ", o_r1);
echo("outer radius large: ", o_r2);
echo("pitch radius small: ", p_r1);
echo("pitch radius large: ", p_r2);
echo("pin_radius: ", pin_radius);
echo("pin_height: ", pin_height);
echo("bar_width: ", bar_width);
echo("bar_length: ", bar_length);
echo("n1: ", n1);
echo("n2: ", n2);
echo("n3: ", n3);