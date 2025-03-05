document.addEventListener("DOMContentLoaded", function () {
    console.log("DOM fully loaded and parsed");
    
    // Locate sidebar, main content, and navbar
    const sidebar = document.querySelector(".bd-sidebar-primary");
    const content = document.querySelector(".bd-main");
    const navbarEnd = document.querySelector(".navbar-item");

    if (!sidebar) {
        console.error("Sidebar not found. Check your theme structure.");
        return;
    }

    if (!navbarEnd) {
        console.error("Navbar not found. Check your theme structure.");
        return;
    }

    // Hide sidebar initially
    sidebar.classList.add("hidden");
    content.classList.add("full-width");

    // Create the toggle button
    const toggleButton = document.createElement("button");
    toggleButton.innerHTML = "☰"; // Hamburger menu icon
    toggleButton.className = "sidebar-toggle-button";
    toggleButton.title = "Toggle Sidebar";

    // Append button to navbar
    navbarEnd.appendChild(toggleButton);
    console.log("Toggle button added to navbar");

    // // Function to check screen size and hide/show the button
    // function checkScreenSize() {
    //     if (window.innerWidth <= 960) {
    //         toggleButton.style.display = "none"; // Hide button in mobile/tablet mode
    //     } else {
    //         toggleButton.style.display = "block"; // Show button in larger screens
    //     }
    // }

    // // Run the function when the page loads
    // checkScreenSize();

    // // Also run the function when the window is resized
    // window.addEventListener("resize", checkScreenSize);

    // Add click event to toggle sidebar visibility
    toggleButton.addEventListener("click", function () {
        console.log("Toggle button clicked");
        sidebar.classList.toggle("hidden");
        content.classList.toggle("full-width");
    });
});
