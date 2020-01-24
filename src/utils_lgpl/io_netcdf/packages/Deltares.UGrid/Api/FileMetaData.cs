﻿using Deltares.UGrid.Entities;
using Deltares.UGrid.Helpers;

namespace Deltares.UGrid.Api
{
    public class FileMetaData
    {
        public FileMetaData(string modelName = "Unknown model", string source = "Unknown Source", string version = "-")
        {
            ModelName = modelName;
            Source = source;
            Version = version;
        }

        /// <summary>
        /// Name of the model for which the file is created
        /// </summary>
        public string ModelName { get; }

        /// <summary>
        /// Source (application) that creates this file
        /// </summary>
        public string Source { get; }

        /// <summary>
        /// Version of the application that creates this file
        /// </summary>
        public string Version { get; }

        internal InteropMetadata CreateMetaData()
        {
            return new InteropMetadata
            {
                institution = "Deltares".ToFixedLengthString(100).ToCharArray(),
                modelname = ModelName.ToFixedLengthString(100).ToCharArray(),
                references = "https://github.com/ugrid-conventions/ugrid-conventions".ToFixedLengthString(100).ToCharArray(),
                source = Source.ToFixedLengthString(100).ToCharArray(),
                version = Version.ToFixedLengthString(100).ToCharArray()
            };
        }
    }
}