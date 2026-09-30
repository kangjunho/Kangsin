(function () {
  'use strict';

  const measurementId = 'G-9S41Y6NWFB';
  const localHosts = new Set(['localhost', '127.0.0.1', '0.0.0.0']);

  if (window.__kangAnalyticsLoaded || location.protocol === 'file:' || localHosts.has(location.hostname)) {
    return;
  }

  window.__kangAnalyticsLoaded = true;
  window.dataLayer = window.dataLayer || [];
  window.gtag = window.gtag || function () {
    window.dataLayer.push(arguments);
  };

  const googleTag = document.createElement('script');
  googleTag.async = true;
  googleTag.src = `https://www.googletagmanager.com/gtag/js?id=${measurementId}`;
  document.head.appendChild(googleTag);

  window.gtag('js', new Date());
  window.gtag('config', measurementId, { send_page_view: false });

  function currentPagePath() {
    return `${location.pathname}${location.search}${location.hash}`;
  }

  window.gtag('event', 'page_view', {
    page_title: document.title,
    page_location: location.href,
    page_path: currentPagePath()
  });

  document.addEventListener('shown.bs.tab', function (event) {
    const target = event.target.getAttribute('data-bs-target') || event.target.getAttribute('href') || '';
    const sectionName = target.replace(/^#/, '') || event.target.textContent.trim();

    window.gtag('event', 'section_view', {
      section_name: sectionName,
      page_title: document.title,
      page_path: currentPagePath()
    });
  });
})();
