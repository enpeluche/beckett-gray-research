const ICONS = {
  PREVIOUS: `<svg class="flip-horizontal" viewBox="0 0 24 24"> <path d="m12 4-1.41 1.41L16.17 11H4v2h12.17l-5.58 5.59L12 20l8-8z"></path> </svg>`,
  NEXT: `<svg viewBox="0 0 24 24"> <path d="m12 4-1.41 1.41L16.17 11H4v2h12.17l-5.58 5.59L12 20l8-8z"></path> </svg>`,
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

export async function loadFooter() {
  const footer = document.getElementById("main-footer");
  if (!footer) return;

  const tocData = await fetchTOC();

  const list = [];
  tocData.forEach((chapter) => {
    list.push({
      title: chapter.title,
      path: `${APP_ROOT}${chapter.path}`,
    });

    if (chapter.children) {
      chapter.children.forEach((subchapter) => {
        list.push({
          title: subchapter.title,
          path: `${APP_ROOT}${subchapter.path}`,
        });
      });
    }
  });

  const pageTitleClean = document.title.toLowerCase().trim();
  const currentIndex = list.findIndex(
    (item) => item.title.toLowerCase().trim() === pageTitleClean
  );

  // Sortir si la page actuelle n'est pas trouvée
  if (currentIndex === -1) {
    console.warn("Titre de page non trouvé dans la liste de navigation.");
    return;
  }

  const previous =
    currentIndex > 0
      ? `
        <div class="previous">
            <a href="${list[currentIndex - 1].path}">
                ${ICONS.PREVIOUS}
                <span class="hhh">
                    ${list[currentIndex - 1].title}
                </span>
            </a>
        </div>`
      : ``;

  const next =
    currentIndex < list.length - 1
      ? `
        <div class="next">
            <a href="${list[currentIndex + 1].path}">
                <span class="hhh">
                    ${list[currentIndex + 1].title}
                </span>
            
                ${ICONS.NEXT}
            </a>
        </div>`
      : ``;

  footer.innerHTML = `${previous} ${next}`;
}
