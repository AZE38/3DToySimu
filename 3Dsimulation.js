import * as THREE from 'https://unpkg.com/three@0.126.1/build/three.module.js';

//import { SpriteText2D, textAlign} from 'https://cdn.jsdelivr.net/npm/three-text2d@0.6.0/lib/index.min.js'
//import { MeshText2D, textAlign } from 'https://cdn.jsdelivr.net/npm/three-text2d@0.6.0/lib/index.min.js'

//import { TextGeometry } from 'https://unpkg.com/three@0.126.1/build/three.module.js'

import { OrbitControls } from 'https://unpkg.com/three@0.126.1/examples/jsm/controls/OrbitControls.js';

//import  Octree  from './octree.js';
//let g = new Octree()

class Vector3 {
    constructor(x, y, z) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
}
class Star {
    constructor(massG, massI, posX, posY, posZ, velX, velY, velZ, color,type) {
        this.type = type;
        this.massG = massG; // Gravitational mass
        this.massI = massI; // Inertial mass
        this.posX = posX; // position x
        this.posY = posY; // position y
        this.posZ = posZ; // position y
        this.velX = velX; // vitesse x
        this.velY = velY; // vitesse y
        this.velZ = velZ; // vitesse y
        this.color = color; // Color of the star
        this.point = null; // point mesh 
    }
      createPoint() {
        let position = new THREE.Vector3(this.posX, this.posY, this.posZ);
        const geometry = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(new Float32Array([position.x, position.y, position.z]), 3));
        const material = new THREE.PointsMaterial({ color: this.color, size: 0.1, sizeAttenuation: true });
        this.point = new THREE.Points(geometry, material);
        this.point.position.set(0, 0, 0);
        scene.add(this.point);
      }
      
      removePoint(){

        scene.remove(this.point);
        this.point.geometry.dispose();
        this.point.material.dispose();
        this.point = null;
      }
      updatePointPositionReplay(starsReplayX,starsReplayY,starsReplayZ){
        this.point.geometry.attributes.position.setXYZ(0,starsReplayX,starsReplayY,starsReplayZ);
        this.point.geometry.attributes.position.needsUpdate = true;
    }
    updatePointPosition(){
        this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        this.point.geometry.attributes.position.needsUpdate = true;
    }
    updatePositionReplay(starsReplayX,starsReplayY,starsReplayZ){ // starsReplay  position at timestep Replay of the particle
        this.posX = starsReplayX;
        this.posY = starsReplayY;
        this.posZ = starsReplayZ;
        this.updatePointPosition()
    }
      

    updatePosition(deltaTime, forceX, forceY, forceZ,BoundR) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;


        this.posX += this.velX * deltaTime;
        this.posY += this.velY * deltaTime;
        this.posZ += this.velZ * deltaTime;


        this.velX += accelX * deltaTime;
        this.velY += accelY * deltaTime;
        this.velZ += accelZ * deltaTime;

        //this.checkBoundaryCondition(speed_init,angleVariance,BoundR);


        this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        this.point.geometry.attributes.position.needsUpdate = true;

    }
    updatePositionLPReplay(deltaTime, forceX, forceY, forceZ,BoundR) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;
        




        this.posX += this.velX * deltaTime + 1/2*accelX*deltaTime*deltaTime ;
        this.posY += this.velY * deltaTime + 1/2*accelY*deltaTime*deltaTime ;
        this.posZ += this.velZ * deltaTime + 1/2*accelZ*deltaTime*deltaTime ;

        this.checkBoundaryCondition(speed_init,angleVariance,BoundR);


        //this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        //this.point.geometry.attributes.position.needsUpdate = true;

    }

    updatePositionLP(deltaTime, forceX, forceY, forceZ,BoundR) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;
        




        this.posX += this.velX * deltaTime + 1/2*accelX*deltaTime*deltaTime ;
        this.posY += this.velY * deltaTime + 1/2*accelY*deltaTime*deltaTime ;
        this.posZ += this.velZ * deltaTime + 1/2*accelZ*deltaTime*deltaTime ;

        this.checkBoundaryCondition(speed_init,angleVariance,BoundR);


        this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        this.point.geometry.attributes.position.needsUpdate = true;

    }
    updateVelocityLP(deltaTime,NewforceX, NewforceY, NewforceZ, forceX, forceY, forceZ) {

        const accelX = forceX / this.massI;
        const accelY = forceY / this.massI;
        const accelZ = forceZ / this.massI;

        const NewaccelX = NewforceX / this.massI;
        const NewaccelY = NewforceY / this.massI;
        const NewaccelZ = NewforceZ / this.massI;
        

        this.velX += 1/2*(accelX+NewaccelX) * deltaTime;
        this.velY += 1/2*(accelY+NewaccelY) * deltaTime;
        this.velZ += 1/2*(accelZ+NewaccelZ) * deltaTime;

    }
    updatePosition2(BoundR) {
        this.checkBoundaryCondition(speed_init,angleVariance,BoundR);
        this.point.geometry.attributes.position.setXYZ(0, this.posX, this.posY, this.posZ);
        this.point.geometry.attributes.position.needsUpdate = true;
       
    }
    checkBoundaryCondition(initialDarkMatterSpeed,angleVariance,BoundRadius) {
        const boundaryRadius = BoundRadius; // Half of the canvas size
        const distanceFromCenter = Math.sqrt(this.posX * this.posX  + this.posY*this.posY + this.posZ*this.posZ);
        

        if (distanceFromCenter > boundaryRadius) {
            //console.log('matter is out',distanceFromCenter,boundaryRadius)

            if (this.type === 'matter') {
                // Move the matter particle outside the canvas to cancel effect far away , -> useless condition , for auto gravitation matter never go outside
                /*this.posX =0
                this.posY =0
                this.posZ =boundaryRadius
                
                this.velX = 0;
                this.velY = 0;
                this.velZ = 0;*/

            } else if (this.type === 'darkMatter') {

                // Place the Nega matter particle at a random position along the circular boundary -> to have effect of "zero divergence" far the galaxy at boundary simulation
                // the size of galaxy need to be small to "dont affect"  Nega matter at boundary , and not have no desirable effet of flow NegaMatter

                // 3D Generate a random angle for theta (0 to pi) and phi (0 to 2*pi)
                /*const theta = Math.random() * Math.PI; // Polar angle
                const phi = Math.random() * 2 * Math.PI; // Azimuthal angle
                const phase = Math.random() * Math.PI;

                // Convert spherical coordinates to Cartesian coordinates for 3D positioning
                this.posX = boundaryRadius* Math.sin(theta) * Math.cos(phi)
                this.posY = boundaryRadius* Math.sin(theta+phase) * Math.sin(phi)
                this.posZ = boundaryRadius* Math.cos(theta+phase)*/


                let posX = 0
                let posY = 0
                let posZ =0

                

                let isgood2 = false;
                let distance=0
                while (isgood2 == false) {
                    posX = Math.random()  * boundaryRadius*2 - boundaryRadius;
                    posY = Math.random()  * boundaryRadius*2 - boundaryRadius;
                    posZ = Math.random()  * boundaryRadius*2 - boundaryRadius;
        
                     distance = Math.sqrt(posX * posX + posY * posY + posZ * posZ);
        
                    if (distance <= boundaryRadius) {
                        isgood2 = true
                    }
                }



                this.posX = ( posX/distance)*boundaryRadius
                this.posY = ( posY/distance)*boundaryRadius
                this.posZ = ( posZ/distance)*boundaryRadius



                   // Calculate the inward direction vector components
                let dx = - this.posX;
                let dy = - this.posY;
                let dz = - this.posZ;

                // Convert vector to spherical coordinates (r, theta, phi)
                let r = Math.sqrt(dx*dx + dy*dy + dz*dz);
                let theta2 = Math.acos(dz / r); // Polar angle from Z-axis
                let phi2 = Math.atan2(dy, dx); // Azimuthal angle from X-axis


                const maxThetaPerturbation = Math.PI / 4; // Adjust this value for less/more variance
                theta2 += (Math.random() - 0.5) * 2 * maxThetaPerturbation; // Randomly perturb theta within the allowed range


                theta2 = Math.max(0, Math.min(Math.PI, theta2));

                // Convert back to Cartesian coordinates
                dx = r * Math.sin(theta2) * Math.cos(phi2);
                dy = r * Math.sin(theta2) * Math.sin(phi2);
                dz = r * Math.cos(theta2);

                // Apply the speed to the direction vector components
                /*this.velX = speed_init * dx;
                this.velY = speed_init * dy;
                this.velZ = speed_init * dz;*/
                let isgood = false;
                let velX = 0;
                let velY = 0;
                let velZ = 0;
                let absVel=0

                while (isgood == false) {
                    velX = Math.random()  * initialDarkMatterSpeed*2 - initialDarkMatterSpeed;
                    velY = Math.random()  * initialDarkMatterSpeed*2 - initialDarkMatterSpeed;
                    velZ = Math.random()  * initialDarkMatterSpeed*2 - initialDarkMatterSpeed;

                    absVel = Math.sqrt(velX * velX + velY * velY + velZ * velZ);

                    if (absVel <= initialDarkMatterSpeed && absVel >= initialDarkMatterSpeed*0.8) {
                        isgood = true
                    }
                }
                this.velX = velX;
                this.velY = velY;
                this.velZ = velZ ;
            }
        }
    }

}
/*function calculateDarkMatterDensity(darkMatterStars, boundaryRadius) {
    let densityNeg = 0;
    darkMatterStars.forEach(star => {
        densityNeg += star.massG;
    });

    const volume =  (4 / 3)*Math.PI * Math.pow(boundaryRadius, 3);
    return densityNeg / volume;
}*/




// Function to generate random values within a specified range
function randomBetween(min, max) {
    return Math.random() * (max - min) + min;
}

function initMatterCluster(numStars,densityNeg,Size,elipseX,elipseZ,speedscale,color='rgba(130, 255, 255, 0.3)') {
    const matterStars = [];
    for (let i = 0; i < numStars; i++) {

        const matterColor = color;

        const mass = densityNeg/numStars; // 10^12/10^4 ( total galaxy solar mass / nmbr star)
        const angle = randomBetween(0, 2 * Math.PI);
        const radius = randomBetween(0.01, Size); // Distance from the center (in pixel)
        const posX =  radius * Math.cos(angle) * elipseX;
        const posY =  radius * Math.sin(angle) ;
        const posZ = randomBetween(0.01, Size/10);


        let ad=speedscale;

        // Velocity for rotational motion
      
        let velX = 0;
        let velY = 0;
        let velZ = 0;

        if ( radius < Size/5 ){
            velX = -Math.sin(angle) * ad*radius;
            velY = Math.cos(angle) * ad*radius;
            velZ = 0
        }
        else{
             velX = -Math.sin(angle) * ad*radius;
             velY = Math.cos(angle) * ad*radius;
             velZ = 0
        }
        
        

        matterStars.push(new Star(mass, Math.abs(mass), posX, posY, posZ, velX, velY,velZ,matterColor,'matter'));
    }
    return matterStars;
}


function initDarkMatterCluster(numStars, speed,densityNeg,size,simulationRadius,color='rgba(255, 0, 0, 1)') {
    const darkMatterStars = [];
    const darkMatterColor =color; 


    for (let i = 0; i < numStars; i++) {
        //const mass = densityNeg/numStars;
        const mass =densityNeg*4/3*Math.PI*simulationRadius*simulationRadius*simulationRadius/numStars

      
        let posX=0 ;
        let posY=0;
        let posZ =0;
        let distance=0;
        let isgood = false;

        while (isgood == false) {
             posX = Math.random()  * simulationRadius*2 - simulationRadius;
             posY = Math.random()  * simulationRadius*2 - simulationRadius;
             posZ = Math.random()  * simulationRadius*2 - simulationRadius;

             distance = Math.sqrt(posX * posX + posY * posY + posZ * posZ);

            if (distance <= simulationRadius && distance >= size) {
                isgood = true
            }
        }


        distance = Math.sqrt(posX * posX + posY * posY + posZ * posZ);

        isgood = false;
        let velX = 0;
        let velY = 0;
        let velZ = 0;
        let absVel=0


     //try to mimic pseudo uniforme random direc
        while (isgood == false) {
            velX = Math.random()  * speed*2 - speed;
            velY = Math.random()  * speed*2 - speed;
            velZ = Math.random()  * speed*2 - speed;

            absVel = Math.sqrt(velX * velX + velY * velY + velZ * velZ);

            if (absVel <= speed && absVel >= speed*0.8) {
                isgood = true
            }
        }
        //console.log('dark position init',radi)
        

        
       // console.log(posX,posY,posZ,velX,velY,velZ)
        darkMatterStars.push(new Star(mass, Math.abs(mass), posX, posY,posZ, velX, velY,velZ,darkMatterColor,'darkMatter'));
    }
    return darkMatterStars;
}
function drawBoundary() {
    const boundaryRadius = canvasWidth/2; // Half of the canvas size
    ctx.beginPath();
    ctx.arc(centerX, centerY, boundaryRadius, 0, 2 * Math.PI);
    ctx.strokeStyle = 'green'; // You can choose any color for the boundary
    ctx.stroke();
}

let computeForcesKernel;
let computeForcesKernel2;

function initializeComputeForcesKernel(numStars) {
    // If a kernel already exists, destroy it
    if (computeForcesKernel) {
        computeForcesKernel.destroy();
    }
    if (computeForcesKernel2) {
        computeForcesKernel2.destroy();
    }
    // Initialize GPU instance if not already initialized or if destroyed
    const gpu = new GPU.GPU();

    // Adjust the output size based on the number of stars
    const outputSize = numStars; //starposi + starnega

    // Create a new kernel with the updated output size
    computeForcesKernel = gpu.createKernel(function(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg) {
        let forceX = 0;
        let forceY = 0;
        let forceZ = 0;
        const myPosX = starPositionsX[this.thread.x];
        const myPosY = starPositionsY[this.thread.x];
        const myPosZ = starPositionsZ[this.thread.x];
        const myMass = starMasses[this.thread.x];

        const f1ratio =0
        //const f2ratio=0
    
        for (let i = 0; i < this.constants.numStars; i++) {
            if (i !== this.thread.x) {
                const dx = starPositionsX[i] - myPosX;
                const dy = starPositionsY[i] - myPosY;
                const dz = starPositionsZ[i] - myPosZ;
                const cutoff = 2; //smouthing gravity short scale epsilon
                let distance= Math.sqrt(dx * dx + dy * dy + dz * dz );

                if(distance >= cutoff) {
                    const forceMagnitude = G * myMass * starMasses[i] / (distance * distance);
                    //f1ratio += forceMagnitude
                    forceX += forceMagnitude * (dx / distance);
                    forceY += forceMagnitude * (dy / distance);
                    forceZ += forceMagnitude * (dz / distance);
                }
            }
        }
    
        // Additional fictive force due to constant Negative matter density outside (homogene density)

        const distanceFromCentere = Math.sqrt(myPosX  * myPosX  + myPosY * myPosY + myPosZ * myPosZ);
        //console.log("distanceFromCenter", distanceFromCentere)
        

        // fictive "shadow of negaMatter outside the simulation" (Big lacune + Nega matter inside = 0 (infinite homogeneous repartition) ) 
        const forceMagnitude = 4/3* Math.PI*G*myMass*densityNeg*distanceFromCentere  // interior of the constante density Shadow NegaBigLacune of simulation



        forceX += forceMagnitude * ((myPosX ) / distanceFromCentere);
        forceY += forceMagnitude * ((myPosY ) / distanceFromCentere);
        forceZ += forceMagnitude * ((myPosZ ) / distanceFromCentere);

        

        //const ratio = f1ratio/f2ratio;
        return [forceX, forceY,forceZ];
    }, {
        constants: { numStars: outputSize },
        output: [outputSize]
    });
}


let isFisrtLoop=true

function updateGalaxyGPU(stars, deltaTime, G, densityNeg, boundaryRadius) {
    //console.log("updateGPU")
    // Prepare data for GPU
    const starMasses = stars.map(star => star.massG);
    const starPositionsX = stars.map(star => star.posX);
    const starPositionsY = stars.map(star => star.posY);
    const starPositionsZ = stars.map(star => star.posZ);

    // First kernel to compute forces

    const forces = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg);
    const forcesX = forces.map(force => force[0]);
    const forcesY = forces.map(force => force[1]);
    const forcesZ = forces.map(force => force[2]);
    
    // Update the stars array with new positions and velocities cpu

    for (let i = 0; i < stars.length; i++) {

       stars[i].updatePosition(deltaTime, forcesX[i], forcesY[i],forcesZ[i],boundaryRadius);

    }
    
}

let forceNow =[]
let forceNew=[]
let forceNowX=[]
let forceNowY=[]
let forceNowZ=[]

function updateGalaxyGPULP(stars, deltaTime, G, densityNeg, boundaryRadius) {
    //console.log("updateGPU")
    // Prepare data for GPU
    for (let i = 0; i < stars.length; i++) {

        stars[i].updatePositionLP(deltaTime, forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);
 
     }




    const starMasses = stars.map(star => star.massG);
    const starPositionsX = stars.map(star => star.posX);
    const starPositionsY = stars.map(star => star.posY);
    const starPositionsZ = stars.map(star => star.posZ);

    // First kernel to compute forces

    const forces = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg);
    forceNew = forces
    const forcesX = forces.map(force => force[0]);
    const forcesY = forces.map(force => force[1]);
    const forcesZ = forces.map(force => force[2]);
    


    
    // Update the stars array with new positions and velocities cpu

    for (let i = 0; i < stars.length; i++) {

       stars[i].updateVelocityLP(deltaTime, forcesX[i], forcesY[i],forcesZ[i],forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);

    }
    
    forceNowX = forcesX;
    forceNowY = forcesY 
    forceNowZ = forcesZ
    
}


let ReplayArray=[]

function updateGalaxyGPULPReplay(stars, deltaTime, G, densityNeg, boundaryRadius) {
    // Prepare data for GPU
    let stepReplayArray=[]
    for (let i = 0; i < stars.length; i++) {

        stars[i].updatePositionLPReplay(deltaTime, forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);

        let arrayStarPos = [stars[i].posX,stars[i].posY,stars[i].posZ]
        stepReplayArray.push(arrayStarPos)
     }
     ReplayArray.push(stepReplayArray)



    const starMasses = stars.map(star => star.massG);
    const starPositionsX = stars.map(star => star.posX);
    const starPositionsY = stars.map(star => star.posY);
    const starPositionsZ = stars.map(star => star.posZ);

    // First kernel to compute forces

    const forces = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg);
    forceNew = forces
    const forcesX = forces.map(force => force[0]);
    const forcesY = forces.map(force => force[1]);
    const forcesZ = forces.map(force => force[2]);
    


    
    // Update the stars array with new positions and velocities cpu

    for (let i = 0; i < stars.length; i++) {

       stars[i].updateVelocityLP(deltaTime, forcesX[i], forcesY[i],forcesZ[i],forceNowX[i], forceNowY[i],forceNowZ[i],boundaryRadius);

    }
    
    forceNowX = forcesX;
    forceNowY = forcesY;
    forceNowZ = forcesZ;
    
}


// Create a scene, camera, and renderer
let scene = new THREE.Scene();
let camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 10000);
let renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Create an OrbitControls object
let controls = new OrbitControls(camera, renderer.domElement);

// Set the position and target of the camera
camera.position.set(0, 0, 800);
controls.target.set(0, 0, 0);
controls.update();

let Radius_SImulation= 1200;
let scale= 0.0005;
let  deltaTime= 1000000;
let  numStarsNeg= 1000;
let   numStarsPos=1000;
let densityNeg= -100000000000;
let TotmassPos= 100000000000;
let hole= 20;
let galactR= 90;
let elips= 1;
let speedRot= 0.00000001;
let speed_init= 0.0000000000000000001;
let colorNeg= '#ff0000';
let colorPos= '#00ff00';
let G= 0;
let timestep = 200;

const angleVariance = Math.PI ;

let stars = [];
let animationFrameId;
let animationFrameId2;
let animationFrameId3;
let animationFrameId4;

let darkMatterDensity=0;

var params = {
    Radius_SImulation: 200,
    scale: 0.0005,
    deltaTime: 0.002, // deltatime simulation in Giga years
    numStarsNeg: 20000,
    numStarsPos:10000,
    densityNeg: -8*Math.pow(10,5), // solar mass / ly^3 (density mass of nega stars)
    TotmassPos: 5*Math.pow(10,11), // ~ milky way solar mass (total mass of the galaxy)
    hole: 25,
    galactR: 25, // 25 unit x 1/scale = 50 000 ly milky way radius
    elips: 1,
    speedRot: 32, //  rad per billion years (10^9 year)
    speed_init: 1200, // ly/Gy
    colorNeg: '#ff0000',
    colorPos: '#00ff00',
    timestep:320,
    initialize: function() {
        Radius_SImulation= params.Radius_SImulation;
        scale = params.scale;
        deltaTime= params.deltaTime;
        numStarsNeg= params.numStarsNeg;
        numStarsPos=params.numStarsPos;
        densityNeg= params.densityNeg;
        TotmassPos= params.TotmassPos;
        hole= params.hole;
        galactR= params.galactR;
        elips= params.elips;
        speedRot= params.speedRot;
        speed_init= params.speed_init;
        colorNeg=params.colorNeg;
        colorPos=params.colorPos;
        G = params.G;
        timestep =params.timestep;

        // Code to initialize the simulation
        //G = 1.56*Math.pow(10,-13) * scale * scale * scale  // Adjusted G from ly^3 / (M☉ * yr^2) to px^3 / (M☉ * yr^2) 

        G = 155800 * scale * scale * scale // unit_engine^3 / (M☉ * Gyr^2) G = 155800 -> ly^3 / (M☉ * Gyr^2) Gyr=10^9 yr
        //console.log("params.G",G)
        let totalstars = numStarsNeg +  numStarsPos
        console.log("G ",G )
        initializeComputeForcesKernel(totalstars);
        initSimulation()
    },
    runRealTime: function() {
        // Code to run the simulation in real time
        //animate2();
        //animateLeap()
        animateLp();
    },
    runWait: function() {
        RunSimuWorker(timestep)
    },
    replay: function() {
        Replay(timestep)
    }
};

var gui = new dat.GUI();


gui.add(params, 'Radius_SImulation', hole+1, 10000).step(10).name('Radius of Simulation in light year').onChange(function(value) {
    Radius_SImulation = parseFloat(value,10);
});

gui.add(params, 'scale', 0, 1).step(0.0001).name('scale 1 unit / 2000 ly = 0.0005').onChange(function(value) {
    scale = parseFloat(value,10);
    
});

gui.add(params, 'deltaTime', 0, 1).step(0.00001).name('deltaTime in Giga-year').onChange(function(value) {
    deltaTime = parseInt(value, 10);
});
gui.add(params, 'timestep', 0, 100000).step(1).name('timestep').onChange(function(value) {
    timestep = parseFloat(value,10);
    
});

gui.add(params, 'numStarsNeg', 0, 10000000).step(1).name('number nega stars').onChange(function(value) {
    numStarsNeg = parseInt(value, 10);
});
gui.add(params, 'numStarsPos', 0, 10000000).step(1).name('number posi stars').onChange(function(value) {
    numStarsPos = parseInt(value, 10);
});

gui.add(params, 'densityNeg', -100000000000, 100000000000).step(100).name('density massNega solar mass / ly^3').onChange(function(value) {
    densityNeg = parseFloat(value,10);
});

gui.add(params, 'TotmassPos',  -10000000000000000, 10000000000000000).step(1000).name('total solar mass Galaxy').onChange(function(value) {
    TotmassPos = parseFloat(value,10);
});

gui.add(params, 'hole', 0, 1000).step(1).name('holei').onChange(function(value) {
    hole = parseInt(value, 10);
});

gui.add(params, 'galactR', 0, Radius_SImulation-1).step(1).name('galactR').onChange(function(value) {
    galactR = parseInt(value, 10);
});

gui.add(params, 'elips', 0, 1).step(0.1).name('elips').onChange(function(value) {
    elips = parseFloat(value,10);
});

gui.add(params, 'speedRot', 0, 1000).step(0.1).name('speedRot rad/Gy').onChange(function(value) {
    speedRot = parseFloat(value,10);
    
    
});

gui.add(params, 'speed_init', 0, 100000).step(1).name('speed_init ly/Gy').onChange(function(value) {
    speed_init = parseFloat(value,10);
    
});

// For colors, you'll need to use a color picker
gui.addColor(params, 'colorNeg').name('colorNeg').onChange(function(value) {
    colorNeg = value;
});

gui.addColor(params, 'colorPos').name('colorPos').onChange(function(value) {
    colorPos = value;
});

gui.add(params, 'initialize');
gui.add(params, 'runRealTime');
gui.add(params, 'runWait');
gui.add(params, 'replay');


function initSimulation() {
    ReplayArray = []
    //initializeComputeForcesKernel(params.numStarsPos+params.numStarsNeg);

    // If there's an existing animation frame, cancel it
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    if (animationFrameId3) {
        cancelAnimationFrame(animationFrameId3);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }
    

    if(stars.length !=0){
        stars.forEach(star => {
        star.removePoint()
        
    })}
    stars = [];
        // Draw the circular boundary
    //drawBoundary();
    // Assume initMatterCluster and initDarkMatterCluster are defined elsewhere
    const matterStars = initMatterCluster(numStarsPos, TotmassPos, galactR, elips,elips, speedRot,colorPos);
    //let speed_init = Math.sqrt(kB * temprature_init);
    const darkMatterStars = initDarkMatterCluster(numStarsNeg, speed_init, densityNeg, hole,Radius_SImulation,colorNeg );
    //darkMatterDensity = calculateDarkMatterDensity(darkMatterStars, Radius_SImulation);

    //densityNeg
    // Combine the stars
    stars = [...matterStars, ...darkMatterStars];
    //console.log("stars.length",stars.length)

    // You can also draw the initial state of stars here if needed
    stars.forEach(star => {
        //star.draw();
        star.createPoint()

        //console.log("star.point",star.point)
    
    });

    const starMasses = stars.map(star => star.massG);
    const starPositionsX = stars.map(star => star.posX);
    const starPositionsY = stars.map(star => star.posY);
    const starPositionsZ = stars.map(star => star.posZ);

    forceNow = computeForcesKernel(starMasses, starPositionsX, starPositionsY,starPositionsZ, G, densityNeg)
    console.log(forceNow)
    forceNowX = forceNow.map(forceNow => forceNow[0]);
    forceNowY = forceNow.map(forceNow => forceNow[1]);
    forceNowZ = forceNow.map(forceNow => forceNow[2]);

    renderer.render(scene, camera);
    animateFree()
}
let forcesTot=[]


function animateLp() {
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }
   
    animationFrameId3=requestAnimationFrame(animateLp);


    updateGalaxyGPULP(stars, deltaTime, G, densityNeg, Radius_SImulation);
    controls.update();
    renderer.render(scene, camera);
    
}
function animate2() {
    animationFrameId2=requestAnimationFrame(animate2);
    updateGalaxyGPU(stars, deltaTime, G, densityNeg, Radius_SImulation);
    controls.update();
    renderer.render(scene, camera);
    
}

function animateFree() {
    if (animationFrameId3) {
        cancelAnimationFrame(animationFrameId3);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }
   
    animationFrameId= requestAnimationFrame(animateFree);
    controls.update();
    renderer.render(scene, camera);

    
}
let ReplaySimu = [];
let stepNow=0



function RunSimuWorker(){
    
    console.log('totalstep',timestep)

    if (animationFrameId3) {
        cancelAnimationFrame(animationFrameId3);
    }
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId4) {
        cancelAnimationFrame(animationFrameId4);
    }

    animationFrameId2= requestAnimationFrame(RunSimuWorker);


    if (stepNow < timestep){
        console.log("stepNow",stepNow)
        updateGalaxyGPULPReplay(stars, deltaTime, G, densityNeg, Radius_SImulation)
    } else{
        console.log("FINISHED ")
    
    }

    stepNow++

    controls.update();
    renderer.render(scene, camera);


}
/*function RunSimuWorker(totaltimestep){
    
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    
    for (let i = 0; i < totaltimestep; i++) {
        updateGalaxyGPULPReplay(stars, deltaTime, G, densityNeg, Radius_SImulation)
    }

}*/


let timeNow=0

function Replay() {
    if (computeForcesKernel) {
        computeForcesKernel.destroy();
    }
    if (computeForcesKernel2) {
        computeForcesKernel2.destroy();
    }
    if (animationFrameId) {
        cancelAnimationFrame(animationFrameId);
    }
    if (animationFrameId2) {
        cancelAnimationFrame(animationFrameId2);
    }
    animationFrameId2= requestAnimationFrame(Replay);

    if (timeNow<timestep){

        for (let j = 0; j < stars.length; j++){
            
            //stars[j].updatePositionReplay(ReplayArray[timeNow][j][0],ReplayArray[timeNow][j][1],ReplayArray[timeNow][j][2])
            stars[j].updatePointPositionReplay(ReplayArray[timeNow][j][0],ReplayArray[timeNow][j][1],ReplayArray[timeNow][j][2])

        }
        timeNow++
    }
    controls.update();
    renderer.render(scene, camera);
}

