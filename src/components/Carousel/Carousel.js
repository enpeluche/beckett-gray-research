export function initCarousel() {
    
    const carousels = document.querySelectorAll(".carousel");
    
    carousels.forEach(carousel => {
        const slides = carousel.querySelectorAll(".carousel img");
        slides[0].classList.add("active");

        carousel.insertAdjacentHTML('beforeend', `
            <button class="nav-btn prev" onclick="changeSlide(-1)">&#10094;</button>
            <button class="nav-btn nextd" onclick="changeSlide(1)">&#10095;</button>
            <div class="counter">
                1 / ${slides.length}
            </div>`);
        
        const previous = carousel.querySelector(".prev");
        const next = carousel.querySelector(".nextd");
        const counter = carousel.querySelector(".counter");
    
        let currentIndex = 0;
    
        const updateSlide = (index) => {
            
            slides.forEach(slide => slide.classList.remove("active"));
    
            if (index >= slides.length) currentIndex = 0;
            else if (index < 0) currentIndex = slides.length - 1;
            else currentIndex = index;
    
            slides[currentIndex].classList.add("active");
    
            if (counter) {
                counter.textContent = `${currentIndex + 1} / ${slides.length}`;
            }
        };
    
        next.addEventListener("click", () => {updateSlide(currentIndex + 1);});
    
        previous.addEventListener("click", () => {updateSlide(currentIndex - 1);});
    });
}
