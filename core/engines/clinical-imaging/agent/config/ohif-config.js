/** OHIF Viewer v3 configuration for Imaging Intelligence Agent.
 *
 * Connects to Orthanc DICOMweb endpoints on port 8042.
 * Served as /app-config.js inside the OHIF container.
 *
 * OHIF is a browser SPA — all DICOMweb requests originate from the
 * user's browser, NOT from the OHIF container.  Therefore the URLs
 * must resolve from the browser's perspective (i.e. the Docker host).
 * We use window.location to build the Orthanc URL dynamically so the
 * config works regardless of whether the host is localhost, an IP, or
 * a DNS name.
 *
 * Author: Adam Jones
 * Date: March 2026
 */

// Build Orthanc base URL from the browser's current hostname so it
// works on localhost, LAN IPs (192.168.x.x), and remote hosts alike.
var defined_orthancBase =
  window.location.protocol + "//" + window.location.hostname + ":8042";

window.config = {
  routerBasename: "/",
  showStudyList: true,
  extensions: [],
  modes: [],
  customizationService: {},
  dataSources: [
    {
      namespace: "@ohif/extension-default.dataSourcesModule.dicomweb",
      sourceName: "orthanc",
      configuration: {
        friendlyName: "Imaging Intelligence Orthanc",
        name: "orthanc",
        wadoUriRoot: defined_orthancBase + "/wado",
        qidoRoot: defined_orthancBase + "/dicom-web",
        wadoRoot: defined_orthancBase + "/dicom-web",
        qidoSupportsIncludeField: false,
        imageRendering: "wadors",
        thumbnailRendering: "wadors",
        enableStudyLazyLoad: true,
        supportsFuzzyMatching: false,
        supportsWildcard: true,
        staticWado: true,
        singlepart: "bulkdata,video",
        bulkDataURI: {
          enabled: true,
          relativeResolution: "studies",
        },
        omitQuotationForMultipartRequest: true,
      },
    },
  ],
  defaultDataSourceName: "orthanc",
};
