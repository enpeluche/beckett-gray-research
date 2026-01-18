const ICONS = {
  OPEN: `<svg viewBox="0 0 24 24"><path d="M19 5v14H5V5zm1.1-2H3.9c-.5 0-.9.4-.9.9v16.2c0 .4.4.9.9.9h16.2c.4 0 .9-.5.9-.9V3.9c0-.5-.5-.9-.9-.9M11 7h6v2h-6zm0 4h6v2h-6zm0 4h6v2h-6zM7 7h2v2H7zm0 4h2v2H7zm0 4h2v2H7z"></path></svg>`,
  CLOSE: `<svg class="CloseIcon" viewBox="0 0 24 24"><path d="M19 6.41 17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z"></path></svg>`,
};

import { APP_ROOT } from "../../main.js";
async function fetchTOC() {
  const paths = [`${APP_ROOT}src/data/toc.json`];

  for (const path of paths) {
    try {
      const res = await fetch(path);
      if (res.ok) return await res.json();
    } catch (e) {}
  }
  return null;
}

export async function loadTableOfContents() {
  // 1. GESTION DE L'ÉTAT (OUVERT/FERMÉ)
  const storedState = localStorage.getItem("is_TOC_open");
  // Par défaut ouvert (true) si pas de sauvegarde
  let is_TOC_open = storedState === null ? true : storedState === "true";

  const TOC = document.getElementById("toc-wrapper");
  if (!TOC) return;

  // CORRECTION IMPORTANTE : On applique l'état visuel tout de suite !
  if (is_TOC_open) {
    TOC.classList.add("is-open");
  }

  const tocData = await fetchTOC();
  if (!tocData) return console.error("Impossible de charger toc.json");

  let chapNum = 0;
  let appendixNum = 0;

  const listHTML = tocData
    .map((chap) => {
      const isActive = document.title === chap.title;
      const isAppendix = chap.type === "appendix";
      const numLabel = isAppendix
        ? chap.isNumbered
          ? `<div class="num">APPENDICE ${++appendixNum}</div>`
          : ``
        : chap.isNumbered
        ? `<div class="num">CHAPITRE ${++chapNum}</div>`
        : ``;
      const activeStyle = isActive
        ? "border:2px solid var(--primary-border-color);"
        : "";

      const subListHTML = (chap.children || [])
        .map((sub, i) => {
          const isSubActive = document.title === sub.title;
          const subStyle = isSubActive ? "background-color: #f7f6f7;" : "";
          return `
                <li style="${subStyle}">
                    <a href="${APP_ROOT}${sub.path}" target="_self" title="${
            sub.title
          }" rel="noopener referrer">
                        <span class="num">${++i}</span>
                        <h4>${sub.title}</h4>
                    </a>
                </li>`;
        })
        .join("");

      return `
            <li>
                <div class="chapter" style="${activeStyle}">
                    <a href="${APP_ROOT}${chap.path}" title="${chap.title}" rel="noopener referrer">
                        <div class="num">${numLabel}</div>
                        <h3>${chap.title}</h3>
                    </a>
                </div>
                <div class="chapter-content">
                    <ul>${subListHTML}</ul>
                </div>
            </li>`;
    })
    .join("");

  // Génération du HTML
  TOC.innerHTML = `        
        <button id='toc-open-button' aria-label='Ouvrir le sommaire'>
            ${ICONS.OPEN}
        </button>
        <div id="toc-body">
            <header class="toc-head">
                <h5>Table des matières</h5>
                <button id="toc-close-button" aria-label="Fermer le sommaire"> 
                    ${ICONS.CLOSE}
                </button>
            </header>
            <ul id="music">
                ${listHTML}
            </ul>
        </div>
    `;

  // 2. FONCTION TOGGLE CORRIGÉE
  const toggle = () => {
    // On bascule la classe
    TOC.classList.toggle("is-open");

    // On vérifie si la classe est présente pour savoir l'état
    const isOpen = TOC.classList.contains("is-open");

    // On sauvegarde
    localStorage.setItem("is_TOC_open", isOpen);
  };

  // Ajout des écouteurs d'événements
  document.getElementById("toc-open-button").onclick = toggle;
  document.getElementById("toc-close-button").onclick = toggle;

  // 3. GESTION DU SCROLL
  const uniqueKey = "scroll_pos";
  const savedScroll = localStorage.getItem(uniqueKey);
  const musicContainer = document.getElementById("music"); // On récupère l'élément une fois

  if (savedScroll && musicContainer) {
    musicContainer.scrollTop = parseInt(savedScroll);
  }

  if (musicContainer) {
    musicContainer.addEventListener("scroll", () => {
      localStorage.setItem(uniqueKey, musicContainer.scrollTop);
    });
  }
}
