$fn=70;
// Parameters:

legwidth = 10;
showinside = 0; // 1: To show the inside structure
enemysize = 2;

thickness = 5;
left_arm_rot_x = 20;
left_arm_rot_y = 0;

kirby(enemysize,legwidth,showinside);

// All the parts

module kirby(enemy_size, leg_width, show_inside){
difference(){
    // Main body
    color([255/256,182/256,193/256,1])sphere(30);
    sphere(25);
    //cylinder(r=10,h=100,center=true);
    if(show_inside){
    cube(50,50,10);
    translate([-49,0,0])cube(50,50,10);
    }
    difference(){
    translate([0,0,-15])cube([leg_width+4*2,70,16],center=true);
    cube([leg_width,70,52],center=true);
    };
    color([229/256,34/256,83/256])translate([30,0,-9])sphere(enemy_size);
     
 }
right_eye();
left_eye();
cheek();
right_arm();
left_arm();
right_leg(leg_width);
left_leg(leg_width);
right_support(leg_width);
left_support(leg_width);
}

module right_leg(leg_width){
translate([0,-20,-26])
color([229/256,34/256,83/256])
 leg(leg_width);
}
 
module leg(leg_width){
union() {
    difference() {
        union() {
            //translate([0, 0, -thickness/2]) cube([leg_width+4*2, 30, thickness], center=true);
            translate([0, 7.5, 10/2]) union() {
                cube([leg_width+4*2, 15, 10], center=true);
                translate([0, 0, 5]) rotate([0, 90, 0]) cylinder(r = 15/2, h = leg_width+4*2, center = true);
            }
        }
        union() {
            translate([0, 7.5, 20/2]) cube([leg_width, 15+2*2-3.9, 20+2*2], center=true);
            //the cube under the concave cylinder
            translate([-leg_width/2-1, 7.5, 10-14/2]) cube([2+0.01, 3.25, 14], center=true);

            //concave cylinder
            translate([-leg_width/2-1, 7.5, 10]) rotate([0, 90, 0]) cylinder(r = 3, h = 2+0.01, center = true);
        }
    }
    translate([0,7.5,-2.5])cube([leg_width+4*2,15,5],center=true);
    //translate([leg_width/2, 7.5, 10]) sphere(r = 3-0.05);
    
    hull(){
    translate([0, -30/2, -5/2]) rotate([0, 90, 0]) cylinder(r = 5/2, h = leg_width+4*2, center = true);
    translate([0,-10,2])sphere(7);
    translate([0, 0, -thickness/2]) cube([leg_width+4*2, 1, thickness], center=true);}
}
}




module left_leg(leg_width){ 
rotate([0,0,180])translate([0,-20,-26])
color([229/256,34/256,83/256])
leg(leg_width);
}


module right_eye(){
translate([25,-6,13])rotate([0,-22,0])eye();
}

module left_eye(){
translate([25,6,13])rotate([0,-22,0])eye();
}
module eye(){
{
    scale([1,1,3.4])color([0,0,0,1])sphere(2.5);
    translate([0.39,0,3.1])rotate([0,-4,0])scale([1,1,2.1])color([1,1,1,1])sphere(2.3);
    translate([0.28,0,-4])rotate([0,-2.5,0])scale([1,1,1.6])color([50/256,92/256,229/256,1])sphere(2);
}
}

module left_arm(){
color([255/256,182/256,193/256,1])translate([0,31,10])rotate([left_arm_rot_x,10,0])hull(){
    translate ([0,10,0])sphere(r=6);
    sphere(r=9);
    }
}

module right_arm(){
    color([255/256,182/256,193/256,1])translate([0,-30,10])hull(){
    translate ([10,-10,10])sphere(r=6);
    sphere(r=9);
    }
}






   
module cheek(){
    color([211/256,88/256,119/256]){
    
    //Left cheek
    translate([25,16,5])scale([1,2,1])sphere(2);
     
    //Right cheek
    mirror([0,-1,0])translate([25,16,5])scale([1,2,1])sphere(2);
    }
 }
    

module right_support(leg_width){
translate([0,7.5-20,-16])rotate([0,-90,0])cylinder(r=3-0.05, h=leg_width+20,center=true);
}

module left_support(leg_width){
// Left support cylinder
translate([0,-7.5+20,-16])rotate([0,-90,0])cylinder(r=3-0.05, h=leg_width+20,center=true);
}


 
 
 
 