// Scroll-triggered fade-in animations
document.addEventListener("DOMContentLoaded", function () {
  var elements = document.querySelectorAll(".fade-in");

  if (!("IntersectionObserver" in window)) {
    // Fallback: show all elements immediately
    elements.forEach(function (el) {
      el.classList.add("visible");
    });
    return;
  }

  var observer = new IntersectionObserver(
    function (entries) {
      entries.forEach(function (entry) {
        if (entry.isIntersecting) {
          entry.target.classList.add("visible");
          observer.unobserve(entry.target);
        }
      });
    },
    { threshold: 0.12 }
  );

  elements.forEach(function (el) {
    observer.observe(el);
  });
});
