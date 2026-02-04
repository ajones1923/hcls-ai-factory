// Scroll-triggered fade-in animations
document.addEventListener("DOMContentLoaded", function () {
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

  document.querySelectorAll(".fade-in").forEach(function (el) {
    observer.observe(el);
  });
});
