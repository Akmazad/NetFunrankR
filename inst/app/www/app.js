// Copy the generated R code to the clipboard.
document.addEventListener("click", function (event) {
  const button = event.target.closest(".nfr-copy");
  if (!button) return;
  const target = document.getElementById(button.dataset.target);
  if (!target || !navigator.clipboard) return;
  navigator.clipboard.writeText(target.innerText).then(function () {
    const original = button.innerHTML;
    button.textContent = "Copied";
    setTimeout(function () { button.innerHTML = original; }, 1500);
  });
});
