function makeReturn() {
    document.querySelectorAll('.return').forEach(elt => elt.prepend(Object.assign(document.createElement('span'), { textContent: "Retourner ", className: "keyword"})));
}

function makeMeta() {
    const strong = document.createElement("strong");

    const algIn = strong.cloneNode(true);
    algIn.textContent = "Entrée :";
    document.querySelectorAll('.in').forEach(elt => elt.prepend(algIn.cloneNode(true)));
    const algOut = document.createElement('strong');
    algOut.textContent = "Sortie :";
    document.querySelectorAll('.out').forEach(elt => elt.prepend(algOut.cloneNode(true)));
}

function makeWhile() {
    const keyword = document.createElement("span");
    keyword.className = "keyword";

    const startWhile = keyword.cloneNode(true);
    startWhile.textContent = "Tant que";
    document.querySelectorAll('.while').forEach(elt => elt.prepend(startWhile.cloneNode(true)));

    // rest'e a faire le mot' cle faire
    const endWhile = keyword.cloneNode(true);
    endWhile.textContent = "Fin tant que";
    endWhile.style = "display: block;";
    document.querySelectorAll('.while').forEach(elt => elt.appendChild(endWhile.cloneNode(true)));
}

function makeFor() {
    const keyword = document.createElement("span");
    keyword.className = "keyword";

    const startFor = keyword.cloneNode(true);
    startFor.textContent = "Pour";
    document.querySelectorAll('.for').forEach(elt => elt.prepend(startFor.cloneNode(true)));
    const endFor = keyword.cloneNode(true);
    endFor.textContent = "Fin pour";
    endFor.style = "display: block;";
    document.querySelectorAll('.for').forEach(elt => elt.appendChild(endFor.cloneNode(true)));
    
}

export async function makeAlgorithm() {
    //for
    makeMeta();

    makeWhile();
    makeFor();
    
    const forHeaderList = document.querySelectorAll('.forHeader');

    const doSpan = document.createElement('span');
    doSpan.classList.add('keyword');
    doSpan.textContent = 'faire';

    const strongA = document.createElement('strong');
    strongA.textContent = ' à';

    forHeaderList.forEach(elt => {

        const forS = elt.querySelector('.forS');
        const forE = elt.querySelector('.forE');

        forS.insertAdjacentElement('afterend', strongA.cloneNode(true));

        forE.insertAdjacentElement('afterend', doSpan.cloneNode(true));

    });

    const keyword = document.createElement("span");
    keyword.className = "keyword";
    const strong = document.createElement("strong");

    
    


    const startIf = keyword.cloneNode(true);
    startIf.textContent = "Si";
    document.querySelectorAll('.if').forEach(elt => elt.prepend(startIf.cloneNode(true)));
    const startElse = keyword.cloneNode(true);
    startElse.textContent = "Sinon";
    document.querySelectorAll('.else').forEach(elt => elt.prepend(startElse.cloneNode(true)));
    
    const endIf = keyword.cloneNode(true);
    endIf.textContent = "Fin si";
    endIf.style = "display: block;";
    document.querySelectorAll('.conditional').forEach(elt => elt.appendChild(endIf.cloneNode(true)));

    makeReturn();
}

