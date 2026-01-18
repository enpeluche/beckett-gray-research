import { APP_ROOT } from "../../main.js";
async function fetchNews() {
  const paths = [`${APP_ROOT}src/data/news.json`];
  for (const path of paths) {
    try {
      const res = await fetch(path);
      if (res.ok) return res.json();
    } catch (e) {}
  }
  return null;
}

export async function createNews() {
  const newsContainer = document.getElementById("news-container");

  if (!newsContainer) return;

  const newsData = await fetchNews();

  const newsList = newsData
    .map((item) => {
      return `
            <li>
                <span class="news-desc">${item.description}</span>
                <time datetime="">${item.date}</time>
            </li>`;
    })
    .join("");

  newsContainer.innerHTML = `
    <h2 style="margin-bottom: 1rem;">Dernières actualités du projet</h2>
        <ul>
            ${newsList}
        </ul>`;
}
