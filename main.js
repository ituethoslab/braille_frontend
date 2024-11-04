import * as THREE from 'three';
import { STLLoader } from 'three/examples/jsm/loaders/STLLoader';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';

function STLViewer(modelPath, elementID) {
    const elem = document.getElementById(elementID);
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(80, elem.clientWidth / elem.clientHeight, 1, 1000);
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(elem.clientWidth, elem.clientHeight);
    elem.appendChild(renderer.domElement);


    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.5;
    controls.enableZoom = true;
    controls.autoRotate = false;
    controls.autoRotateSpeed = 1.5;

    scene.add(new THREE.HemisphereLight(0xffffff, 1.5));

    const loader = new STLLoader();
    loader.load(modelPath, (geometry) => {
        const material = new THREE.MeshPhongMaterial({ color: 0xff8833, specular: 100, shininess: 100 });
        const mesh = new THREE.Mesh(geometry, material);

      //added to attempt to change the view of the object
        //mesh.rotation.x = Math.PI / 2; 
        mesh.rotation.x = Math.PI / -2;
        mesh.rotation.y = 0; 
        mesh.rotation.z = 0; // this is 90 degrees in radians

        scene.add(mesh);

        const middle = new THREE.Vector3();
        geometry.computeBoundingBox();
        geometry.boundingBox.getCenter(middle);
        mesh.geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(-middle.x, -middle.y, -middle.z));

        const largestDimension = Math.max(
            geometry.boundingBox.max.x,
            geometry.boundingBox.max.y,
            geometry.boundingBox.max.z
        );
        camera.position.z = largestDimension * 1.5;

        const animate = () => {
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        };
        animate();
    });
}

// Call STLViewer after the page loads
window.onload = function() {
    STLViewer("test.stl", "scene");
};

