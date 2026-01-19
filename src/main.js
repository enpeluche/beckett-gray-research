import { loadTableOfContents } from "./components/TableOfContents/TableOfContents.js";
import { createNews } from "./components/News/News.js";
import { createProject } from "./components/Project/Project.js";
import { makeArticle } from "./components/Course/Course.js";
import { loadFooter } from "./components/Footer/footer.js";
import { setupKatex } from "./katex-config.js";
import { initCarousel } from "./components/Carousel/Carousel.js";

import { makeAlgorithm } from "./components/Algorithm/Algorithm.js";

document.addEventListener("DOMContentLoaded", () => {
  loadTableOfContents();
  createProject();
  createNews();
  makeArticle();
  setupKatex();
  makeAlgorithm();
  loadFooter();
  initCarousel();
});

// On détecte si on est en local ou sur GitHub Pages
const isProduction = window.location.hostname.includes("github.io");

// Si Prod : on ajoute le nom du repo. Si Local : on reste à la racine.
export const APP_ROOT = isProduction ? "/beckett-gray-research/" : "/";
