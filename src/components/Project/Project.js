import { APP_ROOT } from "../../main.js";

async function fetchProject() {
  const paths = [`${APP_ROOT}src/data/project.json`];

  for (const path of paths) {
    try {
      const res = await fetch(path);
      if (res.ok) return await res.json();
    } catch (e) {
      console.log("ff");
    }
  }
  return null;
}

export async function createProject() {
  const project = document.getElementById("project-container");

  if (!project) return;

  const projectData = await fetchProject();

  let techList = projectData.skills.tech
    .map((tech) => {
      return `
            <span class="tech-tag">
                ${tech}
            </span>
        `;
    })
    .join("");

  let conceptList = projectData.skills.concepts
    .map((concept) => {
      return `
            <span class="concept-tag">
                ${concept}
            </span>
        `;
    })
    .join("");

  let authorsList = projectData.authors
    .map((author) => {
      return `
            <span>
                ${author.name}
            </span>
        `;
    })
    .join("<br>");

  project.insertAdjacentHTML(
    "afterbegin",
    `
        <div class="meta-grid">
            <div>
                <strong class="keyword-item">Contexte</strong>
                <span>${projectData.context}</span>
            </div>
            <div>
                <strong class="keyword-item">Lieu & Université</strong>
                <span>${projectData.location_university}</span>
            </div>
            <div>
                <strong class="keyword-item">Période</strong>
                <span>${projectData.date}</span>
            </div>
        </div>


        <div class="meta-grid">
            <div>
                <strong class="keyword-item">Auteur</strong>
                ${authorsList}
            </div>
            <div>
                <strong class="keyword-item">Encadrant</strong>
                <span>${projectData.supervisor.name}</span>
            </div>
            <div>
                <strong class="keyword-item">Licence</strong>
                <span>
                    <ul style="display:flex; flex-direction: row; align-items: center;">
                        <li>
                            <svg style="width:1rem;" focusable="false" aria-hidden="false" viewBox="0 0 16 16" aria-label="CC"><path d="M8 0.146c-0.068 0-0.132-0.002-0.199 0-4.242 0.108-7.655 3.586-7.655 7.854 0 4.336 3.519 7.854 7.854 7.854s7.854-3.519 7.854-7.854-3.519-7.854-7.854-7.854zM8 14.505c-3.59 0-6.504-2.914-6.504-6.504 0-3.52 2.804-6.394 6.297-6.504 0.070-0.002 0.137 0 0.207 0 3.59 0 6.504 2.914 6.504 6.504s-2.914 6.504-6.504 6.504zM5.331 6.248c0.675 0.015 0.929 0.339 1.026 0.829l1.205 0.004c-0.195-1.128-1.023-1.792-2.201-1.798-1.333-0.006-2.3 1.106-2.304 2.479-0.005 1.466 0.955 2.489 2.289 2.493 0.040 0 0.078 0.002 0.117 0 0.009-0 0.018 0 0.027 0 0.833-0.011 1.966-0.592 1.97-1.727l-1.076-0.003c-0.072 0.491-0.504 0.805-1.003 0.764-0.020-0.001-0.039-0-0.059-0-0.812-0.002-1.131-0.618-1.155-1.524-0.025-0.946 0.525-1.532 1.164-1.517zM10.903 9.289c-0.020-0.001-0.039-0-0.059-0-0.812-0.002-1.131-0.618-1.155-1.524-0.025-0.946 0.525-1.532 1.164-1.517 0.675 0.015 0.929 0.339 1.026 0.829l1.205 0.004c-0.195-1.128-1.023-1.792-2.201-1.798-1.333-0.006-2.3 1.106-2.304 2.479-0.005 1.466 0.955 2.489 2.289 2.493 0.040 0 0.078 0.002 0.117 0 0.009-0 0.018 0 0.027 0 0.833-0.011 1.966-0.592 1.97-1.727l-1.076-0.003c-0.072 0.491-0.505 0.805-1.003 0.764z"></path></svg>
                        </li>
                        <li>
                            <svg style="width:1rem; margin-left:5px;" focusable="false" aria-hidden="false" viewBox="0 0 16 16" aria-label="BY"><path d="M7.798 0.072c-4.282 0.109-7.727 3.62-7.727 7.928 0 4.376 3.552 7.928 7.928 7.928s7.928-3.552 7.928-7.928-3.552-7.928-7.928-7.928c-0.068 0-0.133-0.002-0.201 0zM7.791 1.435c0.070-0.002 0.138 0 0.209 0 3.624 0 6.566 2.941 6.566 6.566s-2.941 6.566-6.566 6.566c-3.624 0-6.566-2.941-6.566-6.566 0-3.553 2.83-6.454 6.356-6.566zM9.146 4.020c0 0.633-0.513 1.146-1.146 1.146s-1.146-0.513-1.146-1.146 0.513-1.146 1.146-1.146c0.633 0 1.146 0.513 1.146 1.146zM6.814 13.126h2.371v-3.81h0.717c0 0-0-2.723 0-3.152-0.001-0.113 0.010-0.338-0.154-0.502s-0.385-0.154-0.502-0.154c-0.117 0-2.376 0-2.493 0s-0.338-0.010-0.502 0.154-0.153 0.389-0.154 0.502c0 0.429 0 3.152 0 3.152h0.717v3.81z"></path></svg>
                        </li>
                        <li>
                            <svg style="width:1rem; margin-left:5px;" focusable="false" aria-hidden="false" viewBox="0 0 16 16" aria-label="SA"><path d="M7.818 0.094c-4.263 0.108-7.692 3.604-7.692 7.893 0 4.357 3.536 7.893 7.893 7.893s7.893-3.536 7.893-7.893c0-4.357-3.536-7.893-7.893-7.893-0.068 0-0.133-0.002-0.2 0zM7.81 1.451c0.070-0.002 0.138 0 0.208 0 3.608 0 6.536 2.928 6.536 6.536s-2.928 6.536-6.536 6.536c-3.608 0-6.536-2.928-6.536-6.536-0-3.537 2.817-6.425 6.328-6.536zM4.594 6.701h-0.676l1.534 1.719 1.534-1.719h-0.553c0.12-0.839 0.857-1.378 1.71-1.31 0.034 0.002 0.066 0 0.1 0 1.388 0 1.935 1.050 1.981 2.597 0.047 1.615-0.889 2.619-1.981 2.597-1.153-0.023-1.589-0.574-1.757-1.41h-2.058c0.339 1.926 1.757 3.056 3.769 3.060 2.278 0.004 3.923-1.902 3.923-4.247 0-2.505-1.645-4.247-3.923-4.247-0.068 0-0.133-0.004-0.2 0-0.016 0.001-0.030-0.001-0.046 0-1.423 0.024-3.355 1.021-3.355 2.96z"></path></svg>
                        </li>
                    </ul>
                </span>
            </div>            
        </div>
    `
  );

  project.insertAdjacentHTML(
    "beforeend",
    `
        <div>
            <strong class="keyword-item">Technologies & Concepts</strong>
            ${techList}
            ${conceptList}
        </div>
    `
  );

  project.insertAdjacentHTML(
    "afterend",
    `
        <div style="display: flex; flex-direction: row;justify-content:center; gap:2rem; align-items: center;">
            <div>
                <a href="${projectData.github}" target="_blank" class="github-link">
                <svg viewBox="0 0 24 24" aria-hidden="true" style="fill:white;margin: 0.1rem;" focusable="false">
                <path
                    d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z" />
            </svg>
                </a>
            </div>
        
            <div>
                <a href="${projectData.intro_folder}" class="lecture">Accéder au cours !</a>
            </div>

            <div>
                <a class="report" href="${projectData.pdf}" target="_blank" download>
                    <img src="/src/img/paper.svg" alt="">
                </a>
            </div>
        </div>
        `
  );
}
