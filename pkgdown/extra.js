// Language tabs: picking R or Python in one tabset switches every R/Python tabset on the
// page, and the choice is remembered for the next page and the next visit.
(function () {
  var KEY = "denmark-language";
  var LANGS = ["R", "Python"];

  function triggers() {
    return Array.prototype.slice.call(
      document.querySelectorAll('[data-bs-toggle="tab"], [data-bs-toggle="pill"]'));
  }
  function label(el) { return (el.textContent || "").trim(); }
  function show(el) {
    if (el.classList.contains("active")) return;
    if (window.bootstrap && window.bootstrap.Tab) {
      window.bootstrap.Tab.getOrCreateInstance(el).show();
    } else {
      el.click();
    }
  }
  function select(lang) {
    triggers().forEach(function (el) { if (label(el) === lang) show(el); });
  }
  function init() {
    var saved = null;
    try { saved = window.localStorage.getItem(KEY); } catch (e) {}
    if (LANGS.indexOf(saved) >= 0) select(saved);
    triggers().forEach(function (el) {
      if (LANGS.indexOf(label(el)) < 0) return;
      el.addEventListener("shown.bs.tab", function () {
        try { window.localStorage.setItem(KEY, label(el)); } catch (e) {}
        select(label(el));
      });
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
