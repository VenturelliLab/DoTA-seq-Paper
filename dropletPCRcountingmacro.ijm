
// Format for macro use
rename("Composite");
run("Duplicate...", "duplicate");
selectWindow("Composite");
run("Split Channels");
selectWindow("C3-Composite");
close();

// Remove noise in C1 for clean threshold
selectWindow("C1-Composite");
rename("C1");
run("Gaussian Blur...", "sigma=4");

// Threshold C1 and convert to binary
setAutoThreshold("Default dark");
setThreshold(1445, 65535); // Adjust lower bound to tune sensitivitey
run("Convert to Mask");
run("Erode"); run("Erode");
run("Dilate"); run("Dilate");

// Remove noise in C2 for clean threshold
selectWindow("C2-Composite");
rename("C2");
run("Gaussian Blur...", "sigma=4");

// Threshold C2 and convert to binary
setAutoThreshold("Default dark");
setThreshold(1145, 65535);	// Adjust lower bound to tune sensitivitey
run("Convert to Mask");
run("Erode"); run("Erode");
run("Dilate"); run("Dilate");

// Find the intersection of C1 and C2
imageCalculator("AND create", "C1","C2");
selectWindow("Result of C1");
rename("C1_AND_C2");

// Count the number of particles in each image
selectWindow("C1");
run("Watershed");
run("Analyze Particles...", "size=100-Infinity show=Nothing summarize");

selectWindow("C2");
run("Watershed");
run("Analyze Particles...", "size=100-Infinity show=Nothing summarize");

selectWindow("C1_AND_C2");
run("Watershed");
run("Analyze Particles...", "size=100-Infinity show=Nothing summarize");

print("Cash me outside");